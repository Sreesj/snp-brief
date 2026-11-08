"""New sources (MyVariant, ClinPGx, SpliceAI, RegulomeDB, constraint, timeline), the local GWAS
catalogue, carrier rates and the ACMG draft."""

from __future__ import annotations

import gzip
import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from snp_brief.interpret import acmg_draft, genotype_rates, one_in
from snp_brief.models import CallSet, GnomadVariant, PopCount
from snp_brief.sources import extra, gwas_local


def wrap(response):
    return {"retrieved": "2026-11-08T00:00:00+00:00", "request": {}, "response": response}


class ParserTests(unittest.TestCase):
    def test_myvariant_scores_and_civic(self) -> None:
        payload = wrap({"hits": [{
            "_id": "chr7:g.140753336A>T",
            "cadd": {"phred": 32},
            "dbnsfp": {"revel": {"score": [0.93, 0.91]}, "alphamissense": {"score": 0.99, "pred": "P"},
                       "sift": {"score": [0.0, 0.01], "pred": "D"}, "polyphen2": {"hdiv": {"score": 0.97, "pred": "D"}}, "genename": "BRAF"},
            "snpeff": {"ann": [{"effect": "missense_variant", "putative_impact": "MODERATE", "genename": "BRAF", "hgvs_p": "p.Val600Glu"}]},
            "civic": {"molecularProfiles": [{"evidenceItems": [
                {"disease": {"name": "Melanoma"}, "therapies": [{"name": "Vemurafenib"}], "evidenceType": "PREDICTIVE",
                 "significance": "SENSITIVITYRESPONSE", "evidenceLevel": "A", "evidenceDirection": "SUPPORTS",
                 "source": {"citationId": "21639808"}, "description": "BRAF V600E predicts response."}]}]},
        }]})
        sc = extra.parse_myvariant(payload)
        self.assertEqual((sc.revel, sc.alphamissense, sc.alphamissense_class, sc.cadd_phred), (0.93, 0.99, "likely pathogenic", 32.0))
        self.assertEqual((sc.sift, sc.sift_pred, sc.polyphen_pred), (0.0, "deleterious", "probably damaging"))
        self.assertEqual((sc.consequence, sc.gene, sc.protein_change), ("missense variant", "BRAF", "p.Val600Glu"))
        self.assertEqual(sc.civic[0].therapies, ["Vemurafenib"])
        self.assertEqual(sc.civic[0].level, "A")
        self.assertEqual(extra.parse_myvariant(wrap({"hits": []})).has_scores, False)

    def test_pharmgkb(self) -> None:
        payload = wrap({"data": [
            {"accessionId": "1", "levelOfEvidence": {"term": "3"}, "relatedChemicals": [{"name": "warfarin"}]},
            {"accessionId": "2", "levelOfEvidence": {"term": "1A"}, "relatedChemicals": [{"name": "clopidogrel"}],
             "location": {"genes": [{"symbol": "CYP2C19"}]}, "types": [{"term": "Efficacy"}],
             "allelePhenotypes": [{"allele": "AA", "phenotype": "reduced response to clopidogrel"}]},
        ]})
        res = extra.parse_pharmgkb(payload)
        first = res.annotations[0]
        self.assertEqual((first.level, first.drugs, first.genes, first.categories), ("1A", ["clopidogrel"], ["CYP2C19"], ["Efficacy"]))
        self.assertIn("reduced response", first.summary)

    def test_spliceai_regulome_constraint_timeline(self) -> None:
        sp = extra.parse_spliceai(wrap({"answers": [{"variant": "chr1-1-A-G", "body": {"scores": [
            {"DS_AG": "0.01", "DS_AL": "0.62", "DS_DG": "0.00", "DS_DL": "0.03", "g_name": "CFTR"}]}}]}))
        self.assertEqual((sp.max_delta, sp.kind, sp.gene), (0.62, "acceptor loss", "CFTR"))
        rg = extra.parse_regulome(wrap({"regulome_score": {"ranking": "2b", "probability": "0.71"}, "features": [{}, {}]}))
        self.assertEqual((rg.ranking, rg.probability, rg.feature_count), ("2b", 0.71, 2))
        self.assertIn("binding", rg.meaning)
        c = extra.parse_constraint({"request": {"gene": "F5"}, "response": {"data": {"gene": {"symbol": "F5", "gnomad_constraint": {"oe_lof_upper": 0.9, "pli": 0.0}}}}})
        self.assertEqual((c.gene, c.loeuf), ("F5", 0.9))
        self.assertIn("tolerant", c.verdict)
        t = extra.parse_timeline(wrap({"hitCount": 12, "years": {"2019": 2, "2024": 10}}))
        self.assertEqual((t.total, t.by_year), (12, {2019: 2, 2024: 10}))


class LocalGwasTests(unittest.TestCase):
    def test_build_and_lookup(self) -> None:
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        src = tmp / "assoc.tsv"
        header = ["SNPS", "MAPPED_TRAIT", "DISEASE/TRAIT", "P-VALUE", "STRONGEST SNP-RISK ALLELE", "OR or BETA", "95% CI (TEXT)",
                  "REPORTED GENE(S)", "MAPPED_GENE", "PUBMEDID", "STUDY ACCESSION", "INITIAL SAMPLE SIZE"]
        rows = [
            ["rs12979860", "hepatitis C virus clearance", "HCV clearance", "1E-20", "rs12979860-C", "2.1", "[1.8-2.5]", "IL28B", "IFNL4", "19684573", "GCST1", "1,000 cases"],
            ["rs1; rs2", "haplotype", "x", "1E-9", "?", "", "", "", "", "1", "GCST2", ""],
        ]
        src.write_text("\n".join("\t".join(r) for r in [header] + rows) + "\n", encoding="utf-8")
        out, count = gwas_local.build(src, tmp)
        self.assertEqual(count, 1)
        import os
        from unittest import mock
        with mock.patch.dict(os.environ, {"SNP_BRIEF_USER_CACHE": str(tmp / "cache")}):
            payload = gwas_local.lookup("rs12979860", tmp)
            empty = gwas_local.lookup("rs1", tmp)
        from snp_brief.providers import parse_gwas
        hits = parse_gwas(payload, "rs12979860").hits
        self.assertEqual((hits[0].trait, hits[0].risk_allele, hits[0].pmid), ("hepatitis C virus clearance", "C", "19684573"))
        self.assertIn("odds ratio 2.1", hits[0].effect)
        self.assertEqual(empty["response"]["ensembl_phenotypes"], [])


def item_with(af_pops, *, revel=None, consequence="missense_variant", splice=None, loeuf=None):
    pops = [PopCount(pid, ac, an) for pid, ac, an in af_pops]
    callset = CallSet("gnomAD v4 exomes", sum(p.ac for p in pops), sum(p.an for p in pops), None, pops)
    variant = GnomadVariant("1-1-A-G", "GRCh38", "1", 1, "A", "G", [callset])
    return SimpleNamespace(
        variants=[variant], calls={},
        scores=extra.PredictionScores(True, revel=revel, consequence=consequence.replace("_", " ")),
        splice=extra.SpliceResult(True, max_delta=splice, kind="donor loss") if splice is not None else extra.SpliceResult(False),
        constraint=extra.Constraint(True, gene="G", loeuf=loeuf),
        vep=SimpleNamespace(most_severe=consequence),
    )


class InterpretTests(unittest.TestCase):
    def test_carrier_rates(self) -> None:
        self.assertEqual(one_in(0.0434), "1 in 23")
        self.assertEqual(one_in(0.000476), "1 in 2,100")
        rates = genotype_rates(item_with([("afr", 1900, 33000), ("nfe", 20, 100000)]))
        overall = rates[0]
        self.assertTrue(overall.group.startswith("All"))
        self.assertEqual(rates[1].group, "African/African American")
        self.assertAlmostEqual(rates[1].carriers, 1 - (1 - 1900 / 33000) ** 2)

    def test_acmg_codes(self) -> None:
        codes = {c.code: c.strength for c in acmg_draft(item_with([("nfe", 0, 100000)], revel=0.95))}
        self.assertEqual(codes, {"PM2": "Supporting", "PP3": "Strong"})
        codes = {c.code for c in acmg_draft(item_with([("afr", 4000, 30000)], revel=0.1))}
        self.assertEqual(codes, {"BA1", "BP4"})
        codes = {c.code for c in acmg_draft(item_with([("nfe", 1, 100000)], consequence="synonymous_variant", splice=0.02))}
        self.assertEqual(codes, {"PM2", "BP7"})
        review = [c for c in acmg_draft(item_with([("nfe", 0, 100000)], consequence="stop_gained", loeuf=0.2)) if c.code == "PVS1"]
        self.assertEqual(review[0].strength, "To review")


if __name__ == "__main__":
    unittest.main()


class CachedPlusLiveTests(unittest.TestCase):
    def test_cached_snp_gets_new_sources_live_into_user_cache(self) -> None:
        import os
        from unittest import mock

        from snp_brief import pipeline

        root = Path(__file__).resolve().parents[1]
        user = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, user)
        self.addCleanup(pipeline._RECENT_FAIL.clear)
        repo = root / "tests" / "cache"
        assert not (repo / "rs6025" / "myvariant.json").exists()
        live = []

        def fake_myvariant(rsid, timeout=20):
            live.append(rsid)
            return {"retrieved": "2026-11-08T00:00:00+00:00", "request": {}, "response": {"hits": [{"cadd": {"phred": 23.4}}]}}

        import contextlib, io
        with mock.patch.object(pipeline.extra, "fetch_myvariant", fake_myvariant), \
             mock.patch.dict(os.environ, {"SNP_BRIEF_EXTRAS": "on"}), contextlib.redirect_stderr(io.StringIO()), \
             mock.patch("snp_brief.http_client.request_json", side_effect=RuntimeError("offline in tests")):
            item = pipeline.load_variant("rs6025", repo_cache=repo, user_cache=user, offline=False, budget=10)
        self.assertEqual(live, ["rs6025"])
        self.assertEqual(item.scores.cadd_phred, 23.4)
        self.assertEqual(item.calls["clinvar"].origin, "cache", "cached sources are not fetched again")
        self.assertTrue((user / "rs6025" / "myvariant.json").is_file(), "the new source is saved to the user cache")
        self.assertFalse((repo / "rs6025" / "myvariant.json").exists(), "the checked-in cache is left alone")


class AcmgConflictTests(unittest.TestCase):
    def test_common_pathogenic_allele_is_flagged(self) -> None:
        import contextlib, io, os
        from unittest import mock

        from snp_brief.pipeline import gather

        cache = Path(__file__).resolve().parents[1] / "tests" / "cache"
        with mock.patch.dict(os.environ, {"SNP_BRIEF_EXTRAS": "on"}), contextlib.redirect_stderr(io.StringIO()):
            item = gather("rs6025", cache_dir=cache, refresh=False, offline=True)
        codes = {c.code: c for c in acmg_draft(item)}
        self.assertIn("BS1", codes)
        self.assertEqual(codes["ClinVar"].strength, "Check first")
        self.assertIn("Pathogenic", codes["ClinVar"].evidence)
