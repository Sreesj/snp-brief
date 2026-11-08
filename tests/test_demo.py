"""Demo page: pasted input, per-field status, populations, sources, and the page itself."""

from __future__ import annotations

import contextlib
import io
import logging
import os
import re
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from streamlit.testing.v1 import AppTest

from snp_brief import cli
from snp_brief.demo.view import (
    FIELDS,
    briefs_text,
    paper_links,
    parse_rsid_text,
    source_rows,
    statements,
    summary_row,
    unresolved_populations,
)
from snp_brief.pipeline import gather

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "tests" / "cache"  # frozen test snapshot, independent of the demo list
APP = ROOT / "src" / "snp_brief" / "demo" / "app.py"
FIXTURE_IDS = (ROOT / "fixture" / "snps.txt").read_text(encoding="utf-8").split()
TEST_IDS = [  # SNPs in the frozen tests/cache snapshot
    "rs6025",
    "rs334",
    "rs1800562",
    "rs1799963",
    "rs4988235",
    "rs1815739",
    "rs1801133",
    "rs1426654",
    "rs12913832",
    "rs6054257",
]


def cached(rsid: str):
    return gather(rsid, cache_dir=CACHE, refresh=False, offline=True)


def sources_by_name(rsid: str) -> dict[str, dict[str, str]]:
    return {row["Source"]: row for row in source_rows(cached(rsid))}


class InputTests(unittest.TestCase):
    def test_fixture_text_gives_the_pinned_ids(self) -> None:
        ids, rejected = parse_rsid_text((ROOT / "fixture" / "snps.txt").read_text(encoding="utf-8"))
        self.assertEqual(ids, FIXTURE_IDS)
        self.assertEqual(rejected, [])

    def test_pasted_list_accepts_separators_and_reports_bad_tokens(self) -> None:
        ids, rejected = parse_rsid_text("# demo list\nrs6025, RS334;rs6025\n\nrs12 x rs0  # trailing note\n")
        self.assertEqual(ids, ["rs6025", "rs334", "rs12"])
        self.assertEqual(rejected, ["x", "rs0"])


class StatusTests(unittest.TestCase):
    def test_statement_status_follows_the_brief_wording(self) -> None:
        found = statements(
            [
                "rs1 is GRCh38 chr1:100 (dbSNP).",
                "Allele T>A: Function: unresolved. missense in A; stop gained in A.",
                "Allele T>G: Function: not reported.",
                "Frequency for G>A: gnomAD v4 genomes allele frequency 0.1 (1/10 alleles).",
                "Disease relevance: no clinical significance was found (ClinVar).",
                "Papers: Smith J, 2020, A title, PMID:1 (Europe PMC).",
            ]
        )
        self.assertEqual(
            [(item.field, item.allele, item.status) for item in found],
            [
                ("Function", "T>A", "unresolved"),
                ("Function", "T>G", "blank"),
                ("Frequency", "G>A", "sourced"),
                ("Disease relevance", "", "blank"),
            ],
        )

    def test_actn3_is_benign_with_ancestry_spread(self) -> None:
        item = cached("rs1815739")
        self.assertEqual(
            summary_row(item),
            {
                "rsID": "rs1815739",
                "Gene": "ACTN3",
                "Function": "sourced",
                "Frequency": "sourced",
                "Disease relevance": "sourced",
            },
        )
        disease = next(sentence for sentence in item.sentences if sentence.startswith("Disease relevance:"))
        self.assertTrue(disease.startswith("Disease relevance: Benign"))
        self.assertIn("Lower-review submissions that disagree: Pathogenic", disease)

    def test_factor_v_leiden_is_pathogenic(self) -> None:
        item = cached("rs6025")
        self.assertEqual(summary_row(item)["Disease relevance"], "sourced")
        disease = next(sentence for sentence in item.sentences if sentence.startswith("Disease relevance:"))
        self.assertTrue(disease.startswith("Disease relevance: Pathogenic"))
        self.assertIn("drug response", disease)
        self.assertIn("Thrombophilia", disease)
        self.assertIn("ClinVar SCV", disease)

    def test_mthfr_stays_unresolved(self) -> None:
        self.assertEqual(summary_row(cached("rs1801133"))["Disease relevance"], "unresolved")

    def test_uniprot_fills_function_and_disease_stays_blank(self) -> None:
        self.assertEqual(
            summary_row(cached("rs6054257")),
            {
                "rsID": "rs6054257",
                "Gene": "DEFB125",
                "Function": "sourced",
                "Frequency": "sourced",
                "Disease relevance": "blank",
            },
        )

    def test_every_fixture_row_agrees_with_its_paragraph(self) -> None:
        for rsid in TEST_IDS:
            item = cached(rsid)
            self.assertEqual(" ".join(item.sentences), item.text)
            row = summary_row(item)
            for name in FIELDS:
                marked = re.search(rf"{name}(?: for [^:]+)?: unresolved\.", item.text) is not None
                self.assertEqual(row[name] == "unresolved", marked, f"{rsid} {name}")

    def test_id_without_cache_is_blank_offline(self) -> None:
        with mock.patch.dict(os.environ, {"SNP_BRIEF_EXTRAS": "on"}), contextlib.redirect_stderr(io.StringIO()):
            item = gather("rs999999999", cache_dir=CACHE, refresh=False, offline=True)
        row = summary_row(item)
        self.assertEqual([row[name] for name in FIELDS], ["blank", "blank", "blank"])
        answered = {entry["Source"]: entry["Answered"] for entry in source_rows(item)}
        self.assertEqual(
            answered,
            {
                "ClinVar": "no",
                "gnomAD": "no",
                "GWAS Catalog": "no",
                "GTEx": "no",
                "dbSNP": "no",
                "Ensembl LD": "no",
                "Ensembl VEP": "no",
                "Open Targets": "no",
                "UniProt (fallback)": "not consulted",
                "NCBI Gene": "not consulted",
                "Europe PMC": "no",
                "dbSNP citations": "no",
                "LitVar2": "no",
            },
        )

class DetailTests(unittest.TestCase):
    def test_populations_are_the_api_counts_and_name_the_brief_groups(self) -> None:
        tables = unresolved_populations(cached("rs1815739"))
        self.assertEqual([table.allele for table in tables], ["C>T"])
        table = tables[0]
        self.assertEqual(table.callsets, ["gnomAD v4 exomes", "gnomAD v4 genomes"])
        self.assertEqual(
            table.named,
            {("African/African American", "gnomAD v4 exomes"), ("Admixed American", "gnomAD v4 exomes")},
        )
        african = next(row for row in table.rows if row["Population"] == "African/African American")
        self.assertEqual(african["gnomAD v4 exomes"], "5239/33480")
        self.assertIn("gnomAD v4 exomes allele frequency 0.445748 (651404/1461372 alleles)", table.totals)

    def test_resolved_allele_has_no_population_table(self) -> None:
        self.assertEqual([table.allele for table in unresolved_populations(cached("rs334"))], ["T>A"])

    def test_uniprot_is_marked_as_the_fallback(self) -> None:
        fallback = sources_by_name("rs6054257")["UniProt (fallback)"]
        self.assertEqual(fallback["Answered"], "yes")
        self.assertIn("Q8N687", fallback["Gave"])
        self.assertIn("function only", fallback["Gave"])
        actn3 = sources_by_name("rs1815739")
        self.assertEqual(actn3["UniProt (fallback)"]["Answered"], "not consulted")
        self.assertEqual(actn3["UniProt (fallback)"]["Detail"], "ClinVar reports clinical significance")
        for name in ("dbSNP", "ClinVar", "gnomAD", "Europe PMC"):
            self.assertEqual(actn3[name]["Answered"], "yes")
            self.assertTrue(actn3[name]["Detail"].startswith("cached, retrieved 2026-"))

    def test_papers_link_to_europe_pmc(self) -> None:
        links = paper_links(cached("rs1815739"))
        self.assertEqual([link.url for link in links][0], "https://europepmc.org/article/MED/21542061")
        self.assertEqual(links[0].citation, "Alfred T et al., 2011, PMID:21542061")

    def test_download_is_what_the_cli_prints(self) -> None:
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            listing = Path(tempfile.mkdtemp()) / "ids.txt"
            listing.write_text("\n".join(TEST_IDS) + "\n", encoding="utf-8")
            code = cli.main(["--offline", "--input", str(listing), "--cache", str(CACHE)])
        self.assertEqual(code, 0)
        self.assertEqual(briefs_text([cached(rsid) for rsid in TEST_IDS]), printed.getvalue())


def _detail_page(cache_dir: str, rsid: str) -> None:
    from pathlib import Path

    from snp_brief.demo.app import render_detail
    from snp_brief.pipeline import gather

    render_detail(gather(rsid, cache_dir=Path(cache_dir), refresh=False, offline=True))


class PageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        os.environ["SNP_BRIEF_EXTRAS"] = "off"  # page tests must not call live providers
        logging.getLogger("streamlit.runtime.scriptrunner_utils.script_run_context").addFilter(
            lambda record: "missing ScriptRunContext" not in record.getMessage()
        )

    def test_landing_is_empty_and_search_takes_many_rsids(self) -> None:
        page = AppTest.from_file(str(APP), default_timeout=60).run()
        self.assertFalse(page.exception)
        self.assertFalse([b for b in page.button if (b.key or "").startswith("pick-")], "nothing is loaded on first visit")
        self.assertFalse(page.toggle(key="offline").value)
        page.toggle(key="offline").set_value(True).run()  # cached SNPs only; no network in tests
        page.text_input(key="search_text").input(", ".join(FIXTURE_IDS)).run()
        self.assertFalse(page.exception)
        picks = [button.key for button in page.button if (button.key or "").startswith("pick-")]
        self.assertEqual(len(picks), len(FIXTURE_IDS))  # one entry per SNP (a merged rsID shows its current ID)
        self.assertEqual(len(page.tabs), 3, "the first SNP opens automatically")
        page.button(key=picks[1]).click().run()
        self.assertFalse(page.exception)

    def test_single_search_opens_the_variant_with_the_button(self) -> None:
        page = AppTest.from_file(str(APP), default_timeout=60).run()
        page.toggle(key="offline").set_value(True).run()
        page.text_input(key="search_text").set_value("rs6025")
        page.button(key="search_go_btn").click().run()
        self.assertFalse(page.exception)
        self.assertEqual(len(page.tabs), 3)
        self.assertTrue(any("rs6025" in block.value for block in page.get("html")))

    def test_detail_shows_brief_populations_and_fallback_source(self) -> None:
        page = AppTest.from_function(_detail_page, args=(str(CACHE), "rs6054257"), default_timeout=30).run()
        self.assertFalse(page.exception)
        blocks = " ".join(block.value for block in page.get("html"))
        self.assertIn("rs6054257", blocks)
        self.assertIn("DEFB125", blocks)
        self.assertIn("UniProt Q8N687", blocks)
        self.assertNotIn("passed the fact check", blocks)
        self.assertIn("35083/41528", page.dataframe[0].value.to_string())
        sources = page.table[0].value
        fallback = sources[sources["Source"] == "UniProt (fallback)"].iloc[0]
        self.assertEqual(fallback["Status"], "Loaded")

    def test_module_tabs_run_on_demand(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=True):
            page = AppTest.from_function(_detail_page, args=(str(CACHE), "rs1801133"), default_timeout=30).run()
            self.assertFalse(page.exception)
            self.assertEqual([tab.label for tab in page.tabs], ["Variant Observatory", "Molecular Ripple", "Hypothesis Forge"])
            page.button(key="run-forge-rs1801133").click().run()
            self.assertFalse(page.exception)
            self.assertTrue(any("Open problem" in block.value for block in page.markdown))
            page.button(key="run-ripple-rs1801133").click().run()
            self.assertFalse(page.exception)
            self.assertTrue(any("Missense coding variant" in block.value for block in page.markdown))


if __name__ == "__main__":
    unittest.main()
