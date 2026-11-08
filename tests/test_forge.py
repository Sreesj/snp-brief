"""Hypothesis Forge, Molecular Ripple, the Nebius client settings, and module caches."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from snp_brief.forge import CONFLICTING, INSUFFICIENT, REPORTED, UNCERTAIN, WELL, forge, rank, triage, validate
from snp_brief.llm import NEBIUS_DEFAULT_MODEL, NEBIUS_URL, LLMClient, parse_json
from snp_brief.pipeline import gather
from snp_brief.ripple import OBSERVED, PREDICTED, UNKNOWN, context_key, ripple
from snp_brief.summary import summarize

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "cache"


def cached(rsid: str):
    return gather(rsid, cache_dir=CACHE, refresh=False, offline=True)


class FakeClient:
    def __init__(self, reply: str, model: str = "fake-model") -> None:
        self.reply = reply
        self.model = model
        self.calls = 0

    def chat(self, system: str, user: str, **kwargs) -> str:
        self.calls += 1
        return self.reply


def forge_reply(facts) -> str:
    gene_fact = next(f.id for f in facts if "intron variant" in f.text)
    trait_fact = next(f.id for f in facts if "LACTASE PERSISTENCE" in f.text)
    return json.dumps(
        {
            "open_questions": ["Which gene does the variant regulate?"],
            "hypotheses": [
                {
                    "title": "Enhancer of a nearby gene",
                    "mechanism": "The allele may change an enhancer that controls a neighbouring gene.",
                    "support": [
                        {"text": "rs4988235 is an intron variant in MCM6", "facts": [gene_fact]},
                        {"text": "It is seen in 99.9% of people", "facts": [gene_fact]},
                    ],
                    "gaps": ["No allele-specific expression data", "Target gene not confirmed"],
                    "discriminating_observation": "Expression differs between alleles in intestinal tissue.",
                    "experiment": "Allele-specific expression analysis in public intestinal RNA-seq",
                    "feasibility": "high",
                },
                {
                    "title": "Tags a nearby causal variant",
                    "mechanism": "The SNP could be linked to the variant that acts.",
                    "support": [{"text": "ClinVar records an association with LACTASE PERSISTENCE", "facts": [trait_fact]}],
                    "gaps": ["No fine-mapping"],
                    "discriminating_observation": "The signal is explained by another variant after fine-mapping.",
                    "experiment": "Statistical fine-mapping of the locus",
                    "feasibility": "medium",
                },
                {
                    "title": "Unsupported idea",
                    "mechanism": "Something else.",
                    "support": [{"text": "BRCA1 is involved", "facts": ["F99"]}],
                    "experiment": "Anything",
                    "feasibility": "low",
                },
            ],
        }
    )


class TriageTests(unittest.TestCase):
    def test_fixture_categories(self) -> None:
        expected = {
            "rs6025": WELL,
            "rs334": WELL,
            "rs1800562": WELL,
            "rs1799963": WELL,
            "rs1815739": WELL,
            "rs1426654": WELL,
            "rs1801133": CONFLICTING,
            "rs4988235": UNCERTAIN,
            "rs12913832": REPORTED,
            "rs6054257": INSUFFICIENT,
        }
        for rsid, category in expected.items():
            found, _ = triage(cached(rsid))
            self.assertEqual(found.category, category, rsid)

    def test_well_characterized_gets_no_hypotheses_and_no_llm_call(self) -> None:
        client = FakeClient("{}")
        result = forge(cached("rs6025"), client=client)
        self.assertEqual(result.hypotheses, [])
        self.assertEqual(client.calls, 0)
        self.assertIn("well characterized", result.note)

    def test_insufficient_evidence_does_not_invent(self) -> None:
        client = FakeClient("{}")
        result = forge(cached("rs6054257"), client=client)
        self.assertEqual(client.calls, 0)
        self.assertFalse(result.triage.wants_hypotheses)
        self.assertIn("too thin", result.note)

    def test_reason_cites_real_facts(self) -> None:
        found, facts = triage(cached("rs1801133"))
        ids = {fact.id for fact in facts}
        self.assertTrue(found.fact_ids)
        self.assertTrue(set(found.fact_ids) <= ids)


class HypothesisTests(unittest.TestCase):
    def test_llm_hypotheses_are_validated_ranked_and_cached(self) -> None:
        item = cached("rs4988235")
        _, facts = triage(item)
        client = FakeClient(forge_reply(facts))
        with tempfile.TemporaryDirectory() as tmp:
            result = forge(item, client=client, cache_dir=Path(tmp))
            self.assertEqual([h.title for h in result.hypotheses], ["Enhancer of a nearby gene", "Tags a nearby causal variant"])
            first = result.hypotheses[0]
            self.assertEqual(len(first.support), 1)
            self.assertIn("(ClinVar)", first.support[0]["rendered"])
            self.assertIn("99.9%", first.rejected_support[0]["text"])
            self.assertEqual(result.ranking[0]["tests"], "Enhancer of a nearby gene")
            self.assertEqual(result.writer, "llm:fake-model")
            self.assertTrue((Path(tmp) / "rs4988235" / "forge.json").is_file())
            again = forge(item, client=client, cache_dir=Path(tmp))
            self.assertEqual(client.calls, 1)
            self.assertIn("cached", again.writer)
            # No key in the environment: the stored result is reused, nothing is called.
            with mock.patch.dict(os.environ, {}, clear=True):
                offline = forge(item, client=None, cache_dir=Path(tmp))
            self.assertEqual(len(offline.hypotheses), 2)

    def test_bad_reply_falls_back_to_stage_1(self) -> None:
        result = forge(cached("rs1801133"), client=FakeClient("not json"))
        self.assertEqual(result.hypotheses, [])
        self.assertIn("failed", result.note)

    def test_rank_prefers_consistent_and_supported_before_easy(self) -> None:
        rows = rank(
            [
                {"title": "Easy tag idea", "experiment": "x", "feasibility": "high", "gaps": ["1"], "support": [{}], "consistency": "contradicted"},
                {"title": "Supported", "experiment": "y", "feasibility": "medium", "gaps": ["1"], "support": [{}, {}], "consistency": "consistent"},
                {"title": "Less support", "experiment": "z", "feasibility": "high", "gaps": ["1", "2"], "support": [{}], "consistency": "consistent"},
            ]
        )
        self.assertEqual([row["tests"] for row in rows], ["Supported", "Less support", "Easy tag idea"])

    def test_invented_rsid_is_scrubbed(self) -> None:
        item = cached("rs1801133")
        _, facts = triage(item)
        fact = next(f.id for f in facts if "missense" in f.text)
        out = validate(
            {"hypotheses": [{"title": "LD", "mechanism": "Linked to rs17367504.", "support": [{"text": "rs1801133 is a missense variant in MTHFR", "facts": [fact]}],
                             "experiment": "Condition on rs17367504", "feasibility": "high", "consistency": "partly"}]},
            facts,
        )
        self.assertNotIn("rs17367504", out["hypotheses"][0]["mechanism"])
        self.assertIn("not named in the retrieved evidence", out["hypotheses"][0]["experiment"])

    def test_literature_mechanism_category(self) -> None:
        found, _ = triage(cached("rs12913832"))
        self.assertEqual(found.category, REPORTED)
        self.assertIn("PMID:22234890", found.reason)

    def test_well_characterized_steps_name_dissent(self) -> None:
        found, _ = triage(cached("rs1815739"))
        self.assertIn("lower-review submissions disagree", found.next_steps[0])
        self.assertNotIn("pathogenicity call is reviewed and agreed", " ".join(found.next_steps))

    def test_validate_drops_hypothesis_without_support(self) -> None:
        out = validate({"hypotheses": [{"title": "X", "support": [{"text": "made up", "facts": []}]}]}, [])
        self.assertEqual(out["hypotheses"], [])

    def test_parse_json_tolerates_fences_and_thinking(self) -> None:
        self.assertEqual(parse_json('<think>hm</think>```json\n{"a": 1}\n```'), {"a": 1})


class RippleTests(unittest.TestCase):
    def test_context_from_consequence(self) -> None:
        self.assertEqual(context_key("missense variant"), "missense")
        self.assertEqual(context_key("nonsense"), "nonsense")
        self.assertEqual(context_key("intron variant"), "intron")
        self.assertEqual(context_key(""), "intergenic")

    def test_links_are_labelled(self) -> None:
        result = ripple(cached("rs4988235"))
        self.assertEqual(result.context_key, "intron")
        labels = {link.label for link in result.links}
        self.assertEqual(labels, {OBSERVED, PREDICTED, UNKNOWN})
        self.assertTrue(any("nearest annotation" in link.text for link in result.links))

    def test_without_clinvar_dbsnp_decides_context_and_gene_source(self) -> None:
        result = ripple(cached("rs6054257"))
        self.assertEqual(result.context_key, "utr")
        self.assertTrue(any(link.source == "dbSNP" and "upstream" in link.text for link in result.links))
        self.assertTrue(any("named by dbSNP" in link.text for link in result.links))
        self.assertFalse(any("named by ClinVar" in link.text for link in result.links))

    def test_missense_priorities(self) -> None:
        result = ripple(cached("rs6025"))
        self.assertEqual(result.context_key, "missense")
        self.assertIn("Protein domain at this position", result.priorities)


class ReadOnlyCacheTests(unittest.TestCase):
    def test_saving_to_a_read_only_place_does_not_crash(self) -> None:
        from snp_brief.module_cache import save_result

        stored = save_result(Path("/proc/not-writable"), "rs1", "forge", "k", {"a": 1})
        self.assertEqual(stored["a"], 1)


class ClientTests(unittest.TestCase):
    def test_nebius_is_the_default_with_its_key(self) -> None:
        with mock.patch.dict(os.environ, {"NEBIUS_API_KEY": "k"}, clear=True):
            client = LLMClient.from_env()
        self.assertEqual(client.base_url, NEBIUS_URL)
        self.assertEqual(client.model, NEBIUS_DEFAULT_MODEL)
        self.assertEqual(client.api_key, "k")
        self.assertEqual(client.provider, "Nebius Token Factory")

    def test_forge_model_override(self) -> None:
        with mock.patch.dict(os.environ, {"NEBIUS_API_KEY": "k", "SNP_BRIEF_FORGE_MODEL": "big"}, clear=True):
            client = LLMClient.from_env()
            self.assertEqual(client.for_forge().model, "big")
            self.assertEqual(client.model, NEBIUS_DEFAULT_MODEL)

    def test_no_settings_means_no_client(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertIsNone(LLMClient.from_env())


class SummaryCacheTests(unittest.TestCase):
    def test_llm_summary_is_cached_and_reused_offline(self) -> None:
        item = cached("rs6025")
        facts = summarize(item).facts
        fn = next(f.id for f in facts if f.text.startswith("rs6025 is a missense"))
        client = FakeClient(f"rs6025 is a missense variant in F5 [{fn}].")
        with tempfile.TemporaryDirectory() as tmp:
            first = summarize(item, "llm", client=client, cache_dir=Path(tmp))
            self.assertTrue(first.writer.startswith("llm:fake-model"))
            with mock.patch.dict(os.environ, {}, clear=True):
                offline = summarize(item, "llm", cache_dir=Path(tmp))
            self.assertTrue(offline.writer.endswith("(cached)"))
            self.assertEqual(client.calls, 1)


if __name__ == "__main__":
    unittest.main()
