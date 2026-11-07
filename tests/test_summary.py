"""Plain-English summary: facts, the checker, the template, and the LLM fallback."""

from __future__ import annotations

import unittest
from pathlib import Path

from snp_brief.llm import split_sentences
from snp_brief.pipeline import gather
from snp_brief.summary import Fact, check, render, summarize

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "cache"
FIXTURE_IDS = (ROOT / "fixture" / "snps.txt").read_text(encoding="utf-8").split()

FACTS = [
    Fact("F1", "rs6025 is a missense variant in F5; NM_000130.4(F5):c.1601G>A (p.Arg534Gln) (ClinVar VCV000000642).", "ClinVar"),
    Fact("F2", "C>T allele frequency in gnomAD v4 exomes is 2.2% (31939/1460848 alleles).", "gnomAD"),
    Fact("F3", "Best-reviewed ClinVar classification: Pathogenic.", "ClinVar"),
]


def cached(rsid: str):
    return gather(rsid, cache_dir=CACHE, refresh=False, offline=True)


class FakeClient:
    model = "fake-model"

    def __init__(self, reply: str) -> None:
        self.reply = reply
        self.prompts: list[str] = []

    def chat(self, system: str, user: str) -> str:
        self.prompts.append(user)
        return self.reply


class CheckerTests(unittest.TestCase):
    def test_supported_sentence_passes(self) -> None:
        kept, rejected = check(["rs6025 is a missense variant in F5 [F1].", "It is seen on 2.2% of alleles [F2]."], FACTS)
        self.assertEqual(len(kept), 2)
        self.assertEqual(rejected, [])

    def test_sentence_without_citation_is_rejected(self) -> None:
        _, rejected = check(["rs6025 causes clotting problems."], FACTS)
        self.assertEqual(rejected[0][1], "cites no fact")

    def test_invented_number_is_rejected(self) -> None:
        _, rejected = check(["It is seen on 5% of alleles [F2]."], FACTS)
        self.assertIn("number 5%", rejected[0][1])

    def test_number_from_an_uncited_fact_is_rejected(self) -> None:
        _, rejected = check(["rs6025 in F5 is seen on 2.2% of alleles [F1]."], FACTS)
        self.assertIn("2.2%", rejected[0][1])

    def test_unknown_gene_is_rejected(self) -> None:
        _, rejected = check(["rs6025 also affects BRCA1 [F1]."], FACTS)
        self.assertIn("BRCA1", rejected[0][1])

    def test_unknown_fact_id_is_rejected(self) -> None:
        _, rejected = check(["rs6025 is pathogenic [F9]."], FACTS)
        self.assertIn("F9", rejected[0][1])

    def test_clinical_advice_is_rejected(self) -> None:
        _, rejected = check(["Carriers should seek a diagnosis [F3]."], FACTS)
        self.assertIn("clinical advice", rejected[0][1])

    def test_render_names_sources(self) -> None:
        text = render(["rs6025 is in F5 [F1].", "Common [F2, F3]."], FACTS)
        self.assertEqual(text, "rs6025 is in F5 (ClinVar). Common (gnomAD, ClinVar).")


class TemplateTests(unittest.TestCase):
    def test_every_fixture_template_passes_its_own_check(self) -> None:
        for rsid in FIXTURE_IDS:
            summary = summarize(cached(rsid))
            self.assertEqual(summary.rejected, [], rsid)
            self.assertGreaterEqual(len(summary.sentences), 2, rsid)
            self.assertNotIn("[F", summary.text, rsid)

    def test_factor_v_leiden_reads_plainly(self) -> None:
        text = summarize(cached("rs6025")).text
        self.assertIn("rs6025 is a missense variant in F5 (p.Arg534Gln) (ClinVar).", text)
        self.assertIn("2.2% of alleles in gnomAD v4 exomes", text)
        self.assertIn("ClinVar's best-supported classification is Pathogenic", text)
        self.assertNotIn("unresolved", text)

    def test_mthfr_says_unresolved(self) -> None:
        self.assertIn("Disease relevance is unresolved", summarize(cached("rs1801133")).text)

    def test_association_is_not_called_pathogenic(self) -> None:
        text = summarize(cached("rs4988235")).text
        self.assertIn("association with LACTASE PERSISTENCE, not as a pathogenicity call", text)

    def test_uniprot_fallback_is_function_only(self) -> None:
        text = summarize(cached("rs6054257")).text
        self.assertIn("(UniProt)", text)
        self.assertIn("ClinVar reports no clinical significance", text)


class LLMTests(unittest.TestCase):
    def test_good_llm_output_is_used_and_sees_only_facts(self) -> None:
        item = cached("rs6025")
        summary_facts = summarize(item).facts
        fn = next(f.id for f in summary_facts if f.text.startswith("rs6025 is a missense"))
        client = FakeClient(f"rs6025 is a missense variant in F5 [{fn}]. It changes p.Arg534Gln [{fn}].")
        summary = summarize(item, "llm", client=client)
        self.assertEqual(summary.writer, "llm:fake-model")
        self.assertEqual(len(summary.sentences), 2)
        self.assertIn("[F1]", client.prompts[0])
        self.assertIn("Facts:", client.prompts[0])

    def test_hallucinating_llm_falls_back_to_template(self) -> None:
        client = FakeClient("rs6025 raises clot risk 80-fold [F1]. It is found in 12% of people [F2].")
        summary = summarize(cached("rs6025"), "llm", client=client)
        self.assertEqual(summary.writer, "template")
        self.assertTrue(any("not in the cited facts" in reason for _, reason in summary.rejected))

    def test_failed_call_falls_back(self) -> None:
        class Broken:
            model = "broken"

            def chat(self, system: str, user: str) -> str:
                raise RuntimeError("connection refused")

        summary = summarize(cached("rs6025"), "llm", client=Broken())
        self.assertEqual(summary.writer, "template")

    def test_split_sentences_after_tags(self) -> None:
        parts = split_sentences("<think>x</think>A is B [F1]. C is D [F2, F3].\nE [F4]")
        self.assertEqual(parts, ["A is B [F1].", "C is D [F2, F3].", "E [F4]."])


if __name__ == "__main__":
    unittest.main()
