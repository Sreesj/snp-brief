"""Abstract findings, evidence graph, steering, report, Ripple cache, ranking, NCBI key, freshness."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from snp_brief import http_client
from snp_brief.forge import forge, rank, triage
from snp_brief.graph import evidence_dot
from snp_brief.literature import Abstract, extract_findings, verify
from snp_brief.pipeline import gather
from snp_brief.providers import parse_vep
from snp_brief.report import build_report, load_notes, save_notes
from snp_brief.ripple import ripple, ripple_cached
from snp_brief.summary import summarize

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "tests" / "cache"  # frozen test snapshot, independent of the demo list
ABSTRACT = Abstract(
    "21193851",
    "Association of the LCT-13910C>T polymorphism with obesity",
    "2011",
    "The LCT-13910C>T (rs4988235) polymorphism was associated with higher BMI in 1,091 Spanish adults (odds ratio 1.40).",
)


def cached(rsid: str):
    return gather(rsid, cache_dir=CACHE, refresh=False, offline=True)


class Recorder:
    model = "fake"

    def __init__(self, reply: str) -> None:
        self.reply = reply
        self.users: list[str] = []

    def chat(self, system, user, **kwargs):
        self.users.append(user)
        return self.reply


class FindingTests(unittest.TestCase):
    def test_verified_quote_and_numbers(self) -> None:
        raw = {
            "findings": [
                {"pmid": "21193851", "claim": "Carriers had higher BMI (odds ratio 1.40).", "quote": "associated with higher BMI in 1,091 Spanish adults (odds ratio 1.40)"},
                {"pmid": "21193851", "claim": "Risk doubled (odds ratio 2.0).", "quote": "associated with higher BMI"},
                {"pmid": "21193851", "claim": "Linked to diabetes.", "quote": "linked to type 2 diabetes"},
                {"pmid": "999", "claim": "x", "quote": "y"},
            ]
        }
        kept, rejected = verify(raw, [ABSTRACT])
        self.assertEqual(len(kept), 1)
        reasons = [r["reason"] for r in rejected]
        self.assertIn("number 2.0 is not in the quote", reasons)
        self.assertIn("quote not found verbatim in the abstract", reasons)
        self.assertIn("PMID not among the abstracts", reasons)

    def test_extract_caches_and_feeds_summary(self) -> None:
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        shutil.copytree(CACHE / "rs4988235", tmp / "rs4988235")
        item = gather("rs4988235", cache_dir=tmp, refresh=False, offline=True)
        reply = json.dumps({"findings": [{"pmid": "21193851", "claim": "The variant was associated with higher BMI.", "quote": "was associated with higher BMI"}]})
        client = Recorder(reply)
        with mock.patch("snp_brief.literature.load_abstracts", return_value=[ABSTRACT]):
            out = extract_findings(item, client, tmp)
            again = extract_findings(item, client, tmp)
        self.assertEqual(len(out["findings"]), 1)
        self.assertEqual(len(client.users), 1, "second call uses the cache")
        self.assertIn("(cached)", again["writer"])
        reloaded = gather("rs4988235", cache_dir=tmp, refresh=False, offline=True)
        self.assertEqual(len(reloaded.findings), 1)
        text = summarize(reloaded).text
        self.assertIn("Papers report: The variant was associated with higher BMI (PMID:21193851)", text)


class GraphAndReportTests(unittest.TestCase):
    def test_graph_styles(self) -> None:
        dot = evidence_dot(cached("rs6025"))
        self.assertTrue(dot.startswith("digraph evidence"))
        self.assertIn("F5", dot)
        self.assertIn("style=solid", dot)

    def test_report_contains_every_module(self) -> None:
        item = cached("rs1801133")
        found = forge(item, use_llm=False)
        report = build_report(item, summarize(item), ripple(item), found, None, "check folate data")
        for heading in ("## Summary", "## Evidence brief", "## Molecular Ripple", "## Hypothesis Forge", "## Sources", "## Notebook"):
            self.assertIn(heading, report)
        self.assertIn("not medical advice", report)

    def test_notes_round_trip(self) -> None:
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        self.assertTrue(save_notes(tmp, "rs1", "hello"))
        self.assertEqual(load_notes(tmp, "rs1"), "hello")
        self.assertEqual(load_notes(tmp, "rs2"), "")


class SteerAndRankTests(unittest.TestCase):
    def test_steer_reaches_the_prompt_and_bypasses_cache(self) -> None:
        item = cached("rs1801133")
        client = Recorder('{"hypotheses": []}')
        forge(item, client=client, steer={"rejected": ["LD tag"], "notes": "computational only"})
        self.assertIn("The user rejected these hypotheses", client.users[0])
        self.assertIn("LD tag", client.users[0])
        self.assertIn("computational only", client.users[0])

    def test_separation_ranks_before_support(self) -> None:
        rows = rank(
            [
                {"title": "A", "experiment": "a", "feasibility": "high", "support": [{}, {}, {}], "separates": [], "consistency": "consistent"},
                {"title": "B", "experiment": "b", "feasibility": "high", "support": [{}], "separates": ["A", "C"], "consistency": "consistent"},
                {"title": "C", "experiment": "c", "feasibility": "high", "support": [{}], "against": [{}], "separates": ["A"], "consistency": "consistent"},
            ]
        )
        self.assertEqual([r["tests"] for r in rows], ["B", "A", "C"])


class RippleAndVepTests(unittest.TestCase):
    def test_ripple_result_file(self) -> None:
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        item = cached("rs6025")
        first = ripple_cached(item, tmp)
        self.assertTrue((tmp / "rs6025" / "ripple.json").is_file())
        second = ripple_cached(item, tmp)
        self.assertEqual(first.links, second.links)

    def test_vep_protein_domain(self) -> None:
        vep = parse_vep({"response": [{"most_severe_consequence": "missense_variant", "transcript_consequences": [
            {"gene_symbol": "F5", "canonical": 1, "consequence_terms": ["missense_variant"], "amino_acids": "R/Q", "protein_start": 534,
             "domains": [{"db": "Pfam", "name": "PF00039"}]}
        ]}]})
        self.assertEqual((vep.amino_acids, vep.protein_position, vep.domains), ("R/Q", "534", ["Pfam PF00039"]))


class NcbiKeyAndFreshnessTests(unittest.TestCase):
    def test_api_key_added_to_eutils_only(self) -> None:
        seen = []

        class Resp:
            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def read(self):
                return b"{}"

        def fake_open(request, timeout=0):
            seen.append(request.full_url)
            return Resp()

        with mock.patch.dict(os.environ, {"NCBI_API_KEY": "abc"}), mock.patch("urllib.request.urlopen", side_effect=fake_open):
            http_client.request_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", params={"db": "snp"})
            http_client.request_json("https://rest.ensembl.org/vep/human/id/rs1", params={"x": "1"})
        self.assertIn("api_key=abc", seen[0])
        self.assertNotIn("api_key", seen[1])

    def test_old_user_cache_is_refetched(self) -> None:
        from snp_brief import pipeline

        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        with mock.patch.dict(os.environ, {"SNP_BRIEF_USER_CACHE": str(tmp), "SNP_BRIEF_MAX_AGE_DAYS": "30"}):
            self.assertTrue(pipeline._stale(tmp, {"retrieved": "2020-01-01T00:00:00+00:00"}))
            self.assertFalse(pipeline._stale(CACHE, {"retrieved": "2020-01-01T00:00:00+00:00"}))


if __name__ == "__main__":
    unittest.main()


class MergedAndIndelTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp)
        payload = (ROOT / "tests" / "data" / "rs80357713_dbsnp.json").read_text(encoding="utf-8")
        for rsid in ("rs80357713", "rs80357783"):
            (self.tmp / rsid).mkdir()
            (self.tmp / rsid / "dbsnp.json").write_text(payload, encoding="utf-8")

    def test_merged_rsid_is_followed_and_reported(self) -> None:
        import contextlib
        import io

        with mock.patch.dict(os.environ, {"SNP_BRIEF_EXTRAS": "off"}), contextlib.redirect_stderr(io.StringIO()):
            item = gather("rs80357713", cache_dir=self.tmp, refresh=False, offline=True)
        self.assertEqual(item.rsid, "rs80357783")
        self.assertEqual(item.merged_from, "rs80357713")
        text = summarize(item).text
        self.assertTrue(text.startswith("rs80357713 is now rs80357783 in dbSNP (the IDs were merged)"))
        from snp_brief.demo.view import summary_row

        self.assertEqual(summary_row(item)["rsID"], "rs80357713 → rs80357783")

    def test_insertion_alleles_are_parsed(self) -> None:
        from snp_brief.sources.dbsnp import parse_dbsnp

        db = parse_dbsnp(json.loads((ROOT / "tests" / "data" / "rs80357713_dbsnp.json").read_text(encoding="utf-8")), "rs80357713")
        self.assertEqual(db.alleles, [("T", "TT")])
        self.assertEqual(db.current_rsid, "rs80357783")

    def test_tiny_frequency_is_readable(self) -> None:
        from snp_brief.sources.gnomad import percent

        self.assertEqual(percent(1, 1460000), "under 0.001%")
        self.assertEqual(percent(3, 39666), "0.0076%")

    def test_condition_name_with_diagnostic_is_not_advice(self) -> None:
        from snp_brief.summary import Fact, check

        facts = [Fact("F1", "Best-reviewed ClinVar classification: Pathogenic, conditions Cystic fibrosis diagnostic test.", "ClinVar")]
        kept, rejected = check(["ClinVar calls it Pathogenic for Cystic fibrosis diagnostic test [F1]."], facts)
        self.assertEqual(rejected, [])
        _, rejected = check(["Carriers should get a diagnosis [F1]."], [Fact("F1", "Pathogenic.", "ClinVar")])
        self.assertTrue(rejected)

    def test_uniprot_evidence_lists_are_stripped(self) -> None:
        from snp_brief.summary import _strip_refs

        self.assertEqual(_strip_refs("Plays a role in DNA repair (PubMed:10500182, PubMed:12887909)."), "Plays a role in DNA repair.")
