"""Ranked providers: parsers, parallel gather with a time budget, fallbacks, and use downstream.

Payloads below are stand-ins shaped like each API's documented or observed JSON.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import shutil
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

from snp_brief import pipeline
from snp_brief.forge import TRAIT_ASSOCIATION, triage
from snp_brief.pipeline import gather
from snp_brief.providers import (
    PROVIDERS,
    parse_gene_summary,
    parse_gtex,
    parse_gwas,
    parse_ld,
    parse_litvar,
    parse_open_targets,
    parse_vep,
)
from snp_brief.ripple import PREDICTED, ripple
from snp_brief.summary import summarize

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "cache"
RSID = "rs6054257"  # checked-in SNP with no ClinVar significance; extras are added in a temp copy

VEP = {
    "source": "Ensembl VEP",
    "retrieved": "2026-11-07T00:00:00+00:00",
    "response": [
        {
            "most_severe_consequence": "upstream_gene_variant",
            "variant_class": "SNV",
            "transcript_consequences": [
                {"gene_symbol": "DEFB125", "consequence_terms": ["upstream_gene_variant"], "canonical": 1, "impact": "MODIFIER", "cadd_phred": 3.1}
            ],
            "regulatory_feature_consequences": [{"biotype": "promoter_flanking_region", "consequence_terms": ["regulatory_region_variant"]}],
            "motif_feature_consequences": [{"transcription_factors": ["CTCF"]}],
            "colocated_variants": [{"id": RSID, "pubmed": [111, 222]}],
        }
    ],
}
GWAS_V1 = {
    "source": "GWAS Catalog",
    "retrieved": "2026-11-07T00:00:00+00:00",
    "response": {
        "_embedded": {
            "associations": [
                {
                    "pvalueMantissa": 3,
                    "pvalueExponent": -9,
                    "orPerCopyNum": 1.12,
                    "loci": [{"strongestRiskAlleles": [{"riskAlleleName": f"{RSID}-A"}], "authorReportedGenes": [{"geneName": "DEFB125"}]}],
                    "efoTraits": [{"trait": "body height"}],
                    "study": {"publicationInfo": {"pubmedId": "30000001"}},
                },
                {"pvalueMantissa": 2, "pvalueExponent": -15, "efoTraits": [{"trait": "eosinophil count"}]},
            ]
        }
    },
}
GWAS_V2 = {
    "response": {
        "_embedded": {
            "associations": [
                {"p_value": 4e-10, "efo_traits": [{"efo_trait": "asthma"}], "snp_effect_allele": [f"{RSID}-A"], "mapped_genes": ["DEFB125"], "pubmed_id": 31000001, "beta": 0.05}
            ]
        }
    }
}
GTEX = {"retrieved": "2026-11-07T00:00:00+00:00", "response": {"data": [
    {"geneSymbol": "DEFB125", "tissueSiteDetailId": "Testis", "pValue": 1e-12, "nes": -0.4},
    {"geneSymbol": "DEFB126", "tissueSiteDetailId": "Whole_Blood", "pValue": 1e-6, "nes": 0.2},
]}}
LD = {"retrieved": "2026-11-07T00:00:00+00:00", "response": [
    {"variation1": RSID, "variation2": "rs1111111", "r2": "0.95"},
    {"variation1": RSID, "variation2": "rs2222222", "r2": "0.81"},
    {"variation1": RSID, "variation2": "rs3333333", "r2": "0.20"},
]}
OT = {"retrieved": "2026-11-07T00:00:00+00:00", "response": [{"variant_id": "20_85729_G_A", "body": {"data": {"variant": {"credibleSets": {"rows": [
    {"study": {"traitFromSource": "body height"}, "l2GPredictions": {"rows": [{"score": 0.71, "target": {"approvedSymbol": "DEFB125"}}]}}
]}}}}}]}
LITVAR = {"retrieved": "2026-11-07T00:00:00+00:00", "response": {"pmids": [333, 444, 555], "pmids_count": 3}}
GENE = {"request": {"gene": "DEFB125"}, "retrieved": "2026-11-07T00:00:00+00:00", "response": {"result": {"uids": ["245938"], "245938": {"name": "DEFB125", "description": "defensin beta 125", "summary": "Defensins form a family of antimicrobial peptides."}}}}
CITATIONS = {"retrieved": "2026-11-07T00:00:00+00:00", "response": {"citations": [666]}}


class ParserTests(unittest.TestCase):
    def test_vep(self) -> None:
        vep = parse_vep(VEP)
        self.assertEqual(vep.canonical_gene, "DEFB125")
        self.assertEqual(vep.most_severe, "upstream_gene_variant")
        self.assertEqual(vep.regulatory, ["promoter flanking region"])
        self.assertEqual(vep.motifs, ["CTCF"])
        self.assertEqual(vep.pubmed, ["111", "222"])
        self.assertEqual(vep.cadd_phred, 3.1)

    def test_vep_error_and_empty(self) -> None:
        self.assertFalse(parse_vep(None).retrieved)
        self.assertEqual(parse_vep({"response": {"error": "not found"}}).genes, [])

    def test_gwas_v1_and_v2_shapes(self) -> None:
        v1 = parse_gwas(GWAS_V1).hits
        self.assertEqual([h.trait for h in v1], ["eosinophil count", "body height"])  # sorted by p
        height = v1[1]
        self.assertEqual((height.p_text, height.risk_allele, height.effect, height.genes, height.pmid), ("3e-9", f"{RSID}-A", "odds ratio 1.12", ["DEFB125"], "30000001"))
        v2 = parse_gwas(GWAS_V2).hits
        self.assertEqual((v2[0].trait, v2[0].pmid, v2[0].effect), ("asthma", "31000001", "beta 0.05"))

    def test_gwas_unexpected_shape_is_empty_not_an_error(self) -> None:
        self.assertEqual(parse_gwas({"response": {"weird": 1}}).hits, [])

    def test_gtex_ld_ot_litvar_gene(self) -> None:
        self.assertEqual([e.gene for e in parse_gtex(GTEX)], ["DEFB125", "DEFB126"])
        self.assertEqual([p.rsid for p in parse_ld(LD, RSID)], ["rs1111111", "rs2222222"])
        self.assertEqual(parse_open_targets(OT)[0].gene, "DEFB125")
        self.assertEqual(parse_litvar(LITVAR), ["333", "444", "555"])
        self.assertEqual(parse_gene_summary(GENE, "DEFB125").name, "defensin beta 125")
        self.assertIsNone(parse_gene_summary(GENE, "OTHER"))

    def test_ranks_are_ordered(self) -> None:
        ranks = [p.rank for p in PROVIDERS]
        self.assertEqual(ranks, sorted(ranks))
        self.assertEqual(PROVIDERS[0].key, "clinvar")


class EndToEndTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp)
        shutil.copytree(CACHE / RSID, self.tmp / RSID)
        for name, payload in {
            "ensembl_vep": VEP,
            "gwas": GWAS_V1,
            "gtex": GTEX,
            "ensembl_ld": LD,
            "open_targets": OT,
            "litvar": LITVAR,
            "ncbi_gene": GENE,
            "dbsnp_citations": CITATIONS,
        }.items():
            (self.tmp / RSID / f"{name}.json").write_text(json.dumps(payload), encoding="utf-8")
        patcher = mock.patch.dict(os.environ, {"SNP_BRIEF_EXTRAS": "on"})
        patcher.start()
        self.addCleanup(patcher.stop)
        self.item = gather(RSID, cache_dir=self.tmp, refresh=False, offline=True)

    def test_all_providers_answer_from_cache(self) -> None:
        for key in ("gwas", "gtex", "ensembl_ld", "ensembl_vep", "open_targets", "litvar", "dbsnp_citations", "ncbi_gene"):
            self.assertEqual(self.item.calls[key].origin, "cache", key)

    def test_summary_uses_new_evidence_with_labels(self) -> None:
        summary = summarize(self.item)
        text = summary.text
        self.assertEqual(summary.rejected, [])
        self.assertIn("Genome-wide association studies link it to eosinophil count, body height", text)
        self.assertIn("statistical association, not a clinical classification (GWAS Catalog)", text)
        self.assertIn("In GTEx it is linked to the expression of DEFB125, DEFB126", text)
        self.assertIn("Open Targets", text)
        self.assertIn("promoter flanking region", text)
        sources = {fact.source for fact in summary.facts}
        self.assertTrue({"GWAS Catalog", "GTEx", "Ensembl LD", "Ensembl VEP", "Open Targets", "LitVar2"} <= sources)

    def test_forge_no_longer_dead_ends(self) -> None:
        found, _ = triage(self.item)
        self.assertEqual(found.category, TRAIT_ASSOCIATION)
        self.assertTrue(found.wants_hypotheses)
        self.assertIn("2 variants inherited with it", found.open_problem)

    def test_ripple_shows_eqtl_ld_and_predictions(self) -> None:
        result = ripple(self.item)
        texts = " | ".join(link.text for link in result.links)
        self.assertIn("expression of DEFB125 in Testis (alternate allele lowers expression)", texts)
        self.assertIn("Inherited together with 2 nearby variants", texts)
        self.assertTrue(any(link.label == PREDICTED and "Locus-to-gene" in link.text for link in result.links))
        self.assertTrue(any("CADD PHRED 3.1" in link.text for link in result.links))

    def test_extras_off_switch(self) -> None:
        with mock.patch.dict(os.environ, {"SNP_BRIEF_EXTRAS": "off"}):
            item = gather(RSID, cache_dir=self.tmp, refresh=False, offline=True)
        self.assertNotIn("gwas", item.calls)
        self.assertFalse(item.gwas.hits)


class BudgetTests(unittest.TestCase):
    def test_slow_provider_is_pending_and_finishes_in_background(self) -> None:
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        shutil.copytree(CACHE / RSID, tmp / RSID)
        release = threading.Event()

        def slow_gwas(rsid, timeout=20):
            release.wait(5)
            return GWAS_V1

        quick = {"retrieved": "x", "response": {}}
        with mock.patch.dict(os.environ, {"SNP_BRIEF_EXTRAS": "on"}), \
            mock.patch.object(pipeline, "fetch_gwas", side_effect=slow_gwas), \
            mock.patch.object(pipeline, "fetch_vep", return_value=VEP), \
            mock.patch.object(pipeline, "fetch_citations", return_value=quick), \
            mock.patch.object(pipeline, "fetch_litvar", return_value=quick), \
            mock.patch.object(pipeline, "fetch_ld", return_value=quick), \
            mock.patch.object(pipeline, "fetch_gtex", return_value=quick), \
            mock.patch.object(pipeline, "fetch_open_targets", return_value=quick), \
            mock.patch.object(pipeline, "fetch_gene_summary", return_value=GENE), \
            contextlib.redirect_stderr(io.StringIO()):
            started = time.monotonic()
            item = gather(RSID, cache_dir=tmp, refresh=False, offline=False, budget=0.5)
            elapsed = time.monotonic() - started
            self.assertLess(elapsed, 3)
            self.assertEqual(item.calls["gwas"].origin, "pending")
            self.assertIn("gwas", item.pending)
            self.assertEqual(item.calls["ensembl_vep"].origin, "live")
            release.set()
            for _ in range(50):
                if (tmp / RSID / "gwas.json").is_file():
                    break
                time.sleep(0.1)
            self.assertTrue((tmp / RSID / "gwas.json").is_file(), "late provider saved for next time")
            again = gather(RSID, cache_dir=tmp, refresh=False, offline=True)
            self.assertEqual(again.calls["gwas"].origin, "cache")


if __name__ == "__main__":
    unittest.main()
