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


class IntegrityTests(unittest.TestCase):
    def test_cached_variants_pass_and_bad_values_fail(self) -> None:
        import contextlib, io, os
        from unittest import mock

        from snp_brief.integrity import check
        from snp_brief.pipeline import gather
        from snp_brief.summary import summarize

        cache = Path(__file__).resolve().parents[1] / "tests" / "cache"
        with mock.patch.dict(os.environ, {"SNP_BRIEF_EXTRAS": "on"}), contextlib.redirect_stderr(io.StringIO()):
            for rsid in ("rs6025", "rs1801133", "rs334"):
                item = gather(rsid, cache_dir=cache, refresh=False, offline=True)
                report = check(item, summarize(item))
                self.assertEqual(report.failures, [], rsid)
                self.assertIn("gnomAD counts match the source", report.passed)
            item = gather("rs6025", cache_dir=cache, refresh=False, offline=True)
        item.variants[0].callsets[0].ac += 1  # a parser bug
        item.scores = extra.PredictionScores(True, revel=1.7)
        item.regulome = extra.RegulomeResult(True, ranking="9z")
        failures = " ".join(check(item).failures)
        self.assertIn("but gnomAD returned", failures)
        self.assertIn("REVEL 1.7 is outside 0 to 1", failures)
        self.assertIn("RegulomeDB rank '9z'", failures)


class PValueUnderflowTests(unittest.TestCase):
    def test_tiny_p_value_text_is_valid(self) -> None:
        from snp_brief.integrity import check
        from snp_brief.providers import GwasHit, GwasResult

        item = item_with([("nfe", 1, 100000)])
        item.rsid, item.db, item.records = "rs1", None, []
        item.gwas = GwasResult(True, [GwasHit("C-reactive protein", 0.0, "1.00e-1420"), GwasHit("x", 2.0, "2")])
        warnings = check(item).warnings
        self.assertEqual(len(warnings), 1)
        self.assertIn("p-value 2", warnings[0])


class PValueOrderTests(unittest.TestCase):
    def test_underflowed_p_values_sort_by_their_text(self) -> None:
        from snp_brief.providers import GwasHit, _p_order

        hits = [GwasHit("a", 0.0, "1.00e-1420"), GwasHit("b", 0.0, "2.00e-9237"), GwasHit("c", 1e-8, "1e-8")]
        self.assertEqual([h.trait for h in sorted(hits, key=_p_order)], ["b", "a", "c"])


class StrandTests(unittest.TestCase):
    def test_reverse_strand_risk_allele_is_flipped_but_palindromes_are_not(self) -> None:
        from snp_brief.providers import GwasHit, orient_risk_alleles

        hits = [GwasHit("Type 2 diabetes", 1e-50, "1e-50", "A"), GwasHit("BMI", 1e-9, "1e-9", "T")]
        orient_risk_alleles(hits, [("C", "T")])
        self.assertEqual([(h.risk_allele, h.reported_allele) for h in hits], [("T", "A"), ("T", "")])
        pal = [GwasHit("x", 1e-9, "1e-9", "T")]
        orient_risk_alleles(pal, [("A", "T")])
        self.assertEqual(pal[0].risk_allele, "T")
        pal = [GwasHit("x", 1e-9, "1e-9", "G")]
        orient_risk_alleles(pal, [("C", "G"), ("C", "T")])
        self.assertEqual((pal[0].risk_allele, pal[0].reported_allele), ("G", ""))

    def test_regulome_boolean_features(self) -> None:
        rg = extra.parse_regulome(wrap({"regulome_score": {"ranking": "1f", "probability": "0.55"},
                                        "features": {"ChIP": True, "QTL": True, "Footprint": False}}))
        self.assertEqual(rg.feature_count, 2)


class BorzoiTests(unittest.TestCase):
    def test_parse_and_concordance(self) -> None:
        from snp_brief.providers import Eqtl

        payload = wrap({"model": "Borzoi, 4 replicates", "alleles": [
            {"ref": "C", "alt": "T", "ref_matches_genome": True,
             "genes": [{"gene": "IFNL4", "mean_lfc": -0.12, "max_abs_lfc": 0.4, "distance": 1200,
                        "tissues": {"liver": -0.3, "Whole Blood": 0.02}}],
             "chromatin": [{"assay": "DNASE", "label": "HepG2", "lfc": -0.5}]},
            {"ref": "C", "alt": "G", "ref_matches_genome": True, "genes": []}]})
        bz = extra.parse_borzoi(payload)
        self.assertEqual(bz.alleles[0].genes[0].tissues["liver"], -0.3)
        eqtls = [Eqtl("IFNL4", "Liver", 1e-9, -0.8, "chr19_39248147_C_T_b38"),
                 Eqtl("IFNL4", "Whole Blood", 1e-5, -0.4, "chr19_39248147_C_T_b38"),
                 Eqtl("NOTHERE", "Liver", 1e-5, 0.4, "")]
        conc = extra.concordance(eqtls, bz)
        self.assertEqual([(c.tissue, c.matched_track, c.agrees, c.clear) for c in conc],
                         [("Liver", "liver", True, True), ("Whole Blood", "Whole Blood", False, True)])
        faint = extra.parse_borzoi(wrap({"alleles": [{"ref": "C", "alt": "T", "genes": [
            {"gene": "IFNL4", "mean_lfc": 0.001, "tissues": {"liver": 0.002}}]}]}))
        self.assertFalse(extra.concordance(eqtls[:1], faint)[0].clear)


class ZygosityTests(unittest.TestCase):
    def test_build_modes_and_copy_meaning(self) -> None:
        import contextlib, io, os
        from unittest import mock

        from snp_brief import zygosity
        from snp_brief.pipeline import gather
        from snp_brief.summary import summarize

        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        hpoa = tmp / "phenotype.hpoa"
        hpoa.write_text(
            "#description: test\n"
            "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration\n"
            "OMIM:603903\tSickle cell anemia\t\tHP:0000007\tOMIM:603903\tIEA\t\t\t\t\tI\tx\n"
            "OMIM:603903\tSickle cell anemia\t\tHP:0001903\tOMIM:603903\tIEA\t\t\t\t\tP\tx\n"
            "OMIM:188050\tThrombophilia due to thrombin defect\t\tHP:0000006\tOMIM:188050\tIEA\t\t\t\t\tI\tx\n",
            encoding="utf-8",
        )
        out, rows = zygosity.build(hpoa, tmp)
        self.assertEqual(rows, 2)
        cache = Path(__file__).resolve().parents[1] / "tests" / "cache"
        zygosity._modes.cache_clear()
        self.addCleanup(zygosity._modes.cache_clear)
        with mock.patch.dict(os.environ, {"SNP_BRIEF_ROOT": str(tmp), "SNP_BRIEF_EXTRAS": "on"}), contextlib.redirect_stderr(io.StringIO()):
            item = gather("rs334", cache_dir=cache, refresh=False, offline=True)
            rows = zygosity.check(item)
            main = next(z for z in rows if z.allele == "T>A")
            sickle = next(c for c in main.conditions if c.modes)
            self.assertEqual(sickle.modes, ["autosomal recessive"])
            self.assertIn("carrier", sickle.one_copy)
            self.assertIn("causes disease", sickle.two_copies)
            self.assertEqual(main.change, "p.Glu7Val")
            text = summarize(item).text
        self.assertIn("For the T>A (p.Glu7Val) allele", text)
        self.assertIn("Hb SS disease is autosomal recessive (one copy: carrier, usually unaffected; two copies: causes disease)", text)
        self.assertIn("inheritance is not stated for", text)
        self.assertIn("T>G (p.Glu7Ala): Likely benign", text)

    def test_flags(self) -> None:
        from snp_brief.zygosity import AlleleZygosity, Condition, _flags

        rec = AlleleZygosity("T>A", "", "Pathogenic", "", "", [Condition("x", [], ["autosomal recessive"], "HPO")], 100, 3)
        self.assertIn("3 people with two copies", _flags(rec)[0])
        dom = AlleleZygosity("C>T", "", "Pathogenic", "", "", [Condition("y", [], ["autosomal dominant"], "HPO")], 900, 40)
        self.assertIn("dominant condition", _flags(dom)[0])
        unknown = AlleleZygosity("C>T", "", "Pathogenic", "", "", [Condition("z", [], [], "")], 10, 0)
        self.assertIn("not stated", _flags(unknown)[0])
