"""Conflict rules, allele separation, and the pinned fixture list."""

from __future__ import annotations

import unittest
from pathlib import Path

from snp_brief.brief import assemble_brief
from snp_brief.models import (
    CallSet,
    Classification,
    ClinVarRecord,
    DbSnpRecord,
    GnomadVariant,
    Paper,
    PopCount,
    UniProtHit,
)
from snp_brief.sources.clinvar import parse_clinvar
from snp_brief.sources.dbsnp import parse_dbsnp
from snp_brief.sources.gnomad import disagrees, frequency_phrase

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_IDS = [
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


def db(**kwargs) -> DbSnpRecord:
    base = dict(
        retrieved=True,
        rsid="rs1",
        assembly="GRCh38",
        accession="NC_000001.11",
        chrom="1",
        pos=100,
        alleles=[("A", "T")],
        genes=["GENE"],
        coding_by_alt={"T": "c.1A>T"},
    )
    base.update(kwargs)
    return DbSnpRecord(**base)


def classification(description: str, traits: list[str] | None = None, review: str = "criteria provided") -> Classification:
    return Classification(description, review, traits or [])


def record(**kwargs) -> ClinVarRecord:
    base = dict(
        accession="VCV1",
        title="NM_1.2(GENE):c.1A>T (p.Ala1Thr)",
        genes=["GENE"],
        consequences=["missense variant"],
        cdna="c.1A>T",
        classifications=[classification("Pathogenic", ["Example condition"])],
        no_change=False,
    )
    base.update(kwargs)
    return ClinVarRecord(**base)


def variant(ref: str = "A", alt: str = "T", callsets: list[CallSet] | None = None) -> GnomadVariant:
    if callsets is None:
        callsets = [CallSet("gnomAD v4 genomes", ac=2, an=1000, af=0.002, populations=[])]
    return GnomadVariant("1-100-A-T", "GRCh38", "1", 100, ref, alt, callsets)


class FixtureTests(unittest.TestCase):
    def test_pinned_lists_stay_in_place(self) -> None:
        snps = (ROOT / "fixture" / "snps.txt").read_text(encoding="utf-8")
        self.assertEqual([line.strip() for line in snps.splitlines() if line.strip()], FIXTURE_IDS)
        held = [
            line.strip()
            for line in (ROOT / "fixture" / "held-out.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        self.assertEqual(len(held), 1)
        self.assertNotIn(held[0], FIXTURE_IDS)
        self.assertTrue(held[0].startswith("rs"))


class SourceParserTests(unittest.TestCase):
    def test_dbsnp_maps_genomic_alts_to_coding_changes(self) -> None:
        docsum = (
            "HGVS=NC_000011.10:g.5227002T>A,NC_000011.10:g.5227002T>C,NC_000011.10:g.5227002T>G,"
            "NM_000518.5:c.20A>T,NM_000518.5:c.20A>G,NM_000518.5:c.20A>C|SEQ=[T/A]"
        )
        payload = {
            "response": {
                "result": {
                    "uids": ["334"],
                    "334": {
                        "acc": "NC_000011.10",
                        "chrpos": "11:5227002",
                        "spdi": "NC_000011.10:5227001:T:A,NC_000011.10:5227001:T:C,NC_000011.10:5227001:T:G",
                        "genes": [{"name": "HBB"}],
                        "docsum": docsum,
                    },
                }
            }
        }
        parsed = parse_dbsnp(payload, "rs334")
        self.assertEqual(parsed.assembly, "GRCh38")
        self.assertEqual(parsed.chrom, "11")
        self.assertEqual(parsed.pos, 5227002)
        self.assertEqual(parsed.coding_by_alt, {"A": "c.20A>T", "C": "c.20A>G", "G": "c.20A>C"})

    def test_title_gene_is_used_instead_of_locus_symbols(self) -> None:
        summary = _summary(
            "VCV000015333",
            "334",
            "NM_000518.5(HBB):c.20A>T (p.Glu7Val)",
            "c.20A>T",
            "Pathogenic",
            ["missense variant"],
        )
        summary["gene_sort"] = "LOC1"
        summary["genes"] = [{"symbol": "HBB"}, {"symbol": "LOC106099062"}]
        payload = {"response": {"esummary": {"result": {"uids": ["1"], "1": summary}}}}
        _, records = parse_clinvar(payload, "rs334")
        self.assertEqual(records[0].genes, ["HBB"])

    def test_gene_sort_beats_the_transcript_in_the_title(self) -> None:
        summary = _summary(
            "VCV000007685",
            "4988235",
            "NM_002299.2(LCT):c.-13907C>T",
            "NM_002299.2(LCT):c.-13907C>T",
            "association",
            ["intron variant"],
        )
        summary["gene_sort"] = "MCM6"
        summary["genes"] = [{"symbol": "MCM6"}]
        payload = {"response": {"esummary": {"result": {"uids": ["1"], "1": summary}}}}
        _, records = parse_clinvar(payload, "rs4988235")
        self.assertEqual(records[0].genes, ["MCM6"])
        self.assertEqual(records[0].cdna, "c.-13907C>T")

    def test_clinvar_keeps_this_rsid_and_drops_haplotypes(self) -> None:
        payload = {
            "response": {
                "esummary": {
                    "result": {
                        "uids": ["1", "2", "3"],
                        "1": _summary("VCVKEEP", "10", "NM_1.1(GENE):c.2A>T", "c.2A>T", "Pathogenic", ["missense variant"]),
                        "2": _summary(
                            "VCVHAP",
                            "10",
                            "NM_1.1(GENE):c.[2A>T;9G>A]",
                            "c.[2A>T;9G>A]",
                            "Pathogenic",
                            ["missense variant"],
                        ),
                        "3": _summary("VCVOTHER", "999", "NM_9.1(OTHER):c.3G>A", "c.3G>A", "Benign", ["missense variant"]),
                    }
                }
            }
        }
        retrieved, records = parse_clinvar(payload, "rs10")
        self.assertTrue(retrieved)
        self.assertEqual([item.accession for item in records], ["VCVKEEP"])


class ConflictTests(unittest.TestCase):
    def test_frequency_rule_does_not_invent_a_cutoff_number(self) -> None:
        self.assertFalse(disagrees(0.0218633, 0.0172993))
        self.assertFalse(disagrees(0.40, 0.41))
        self.assertTrue(disagrees(0.004, 0.026))
        self.assertTrue(disagrees(0.0, 0.02))
        self.assertFalse(disagrees(0.0, 0.0005))
        self.assertTrue(disagrees(0.20, 0.35))

    def test_population_disagreement_keeps_both_names_and_counts(self) -> None:
        unresolved, text = frequency_phrase(
            [
                CallSet(
                    "gnomAD v4 genomes",
                    ac=201,
                    an=2000,
                    af=0.1005,
                    populations=[
                        PopCount("afr", 1, 1000),
                        PopCount("nfe", 200, 1000),
                        PopCount("eas", 0, 40),
                    ],
                )
            ]
        )
        self.assertTrue(unresolved)
        self.assertIn("African/African American 1/1000 alleles", text)
        self.assertIn("European (non-Finnish) 200/1000 alleles", text)
        self.assertNotIn("0/40", text)

    def test_close_populations_stay_one_claim(self) -> None:
        unresolved, text = frequency_phrase(
            [CallSet("gnomAD v4 exomes", ac=250, an=2000, af=0.125, populations=[
                PopCount("afr", 100, 1000),
                PopCount("nfe", 150, 1000),
            ])]
        )
        self.assertFalse(unresolved)
        self.assertIn("allele frequency 0.125", text)
        self.assertNotIn("African/African American", text)

    def test_missing_frequency_is_blank_not_zero(self) -> None:
        text = assemble_brief(db(), [record()], [], [], None, clinvar_retrieved=True)
        self.assertIn("Frequency for A>T: not reported.", text)
        self.assertNotIn("allele frequency 0", text)

    def test_one_clinvar_record_stays_on_the_gnomad_allele(self) -> None:
        parsed = db(alleles=[("G", "A"), ("G", "C")], coding_by_alt={}, genes=["MCM6"])
        records = [
            record(
                accession="VCVLCT",
                title="NM_002299.2(LCT):c.-13907C>T",
                genes=["MCM6"],
                consequences=["intron variant"],
                cdna="c.-13907C>T",
                classifications=[classification("association", ["LACTASE PERSISTENCE"])],
            )
        ]
        text = assemble_brief(parsed, records, [variant(ref="G", alt="A")], [], None, clinvar_retrieved=True)
        self.assertEqual(text.count("Disease relevance"), 1)
        self.assertIn("Frequency for G>A", text)
        self.assertIn("association", text)
        self.assertIn("LACTASE PERSISTENCE", text)
        self.assertNotIn("Allele NM_", text)

    def test_wider_population_gap_is_the_one_named(self) -> None:
        unresolved, text = frequency_phrase(
            [
                CallSet(
                    "gnomAD v4 exomes",
                    ac=100,
                    an=4000,
                    af=0.025,
                    populations=[PopCount("afr", 20, 2000), PopCount("nfe", 30, 2000)],
                ),
                CallSet(
                    "gnomAD v4 genomes",
                    ac=50,
                    an=2000,
                    af=0.025,
                    populations=[PopCount("afr", 1, 1000), PopCount("nfe", 400, 1000)],
                ),
            ]
        )
        self.assertTrue(unresolved)
        self.assertIn("gnomAD v4 genomes African/African American 1/1000 alleles", text)
        self.assertIn("gnomAD v4 genomes European (non-Finnish) 400/1000 alleles", text)
        self.assertNotIn("exomes African/African American", text)

    def test_different_alleles_are_not_one_conflict(self) -> None:
        parsed = db(
            alleles=[("T", "A"), ("T", "G")],
            coding_by_alt={"A": "c.20A>T", "G": "c.20A>C"},
            genes=["HBB"],
        )
        records = [
            record(
                accession="VCVSICKLE",
                title="NM_000518.5(HBB):c.20A>T (p.Glu7Val)",
                genes=["HBB"],
                cdna="c.20A>T",
                classifications=[classification("Pathogenic", ["Sickle cell disease"])],
            ),
            record(
                accession="VCVOTHER",
                title="NM_000518.4(HBB):c.20A>C (p.Glu7Ala)",
                genes=["HBB"],
                cdna="c.20A>C",
                classifications=[classification("Likely benign", ["not specified"])],
            ),
        ]
        text = assemble_brief(parsed, records, [], [], None, clinvar_retrieved=True)
        self.assertIn("Allele T>A:", text)
        self.assertIn("Allele T>G:", text)
        self.assertIn("Pathogenic", text)
        self.assertIn("Likely benign", text)
        self.assertNotIn("Disease relevance: unresolved", text)
        self.assertNotIn("not specified", text)

    def test_same_allele_conflict_keeps_both(self) -> None:
        records = [
            record(accession="VCVPATH", classifications=[classification("Pathogenic", ["Condition A"])]),
            record(accession="VCVBENIGN", classifications=[classification("Benign", ["Condition B"])]),
        ]
        text = assemble_brief(db(), records, [], [], None, clinvar_retrieved=True)
        self.assertIn("Disease relevance: unresolved.", text)
        self.assertIn("Pathogenic", text)
        self.assertIn("Benign", text)
        self.assertIn("VCVPATH", text)
        self.assertIn("VCVBENIGN", text)

    def test_conflict_phrase_is_unresolved(self) -> None:
        records = [
            record(classifications=[classification("Conflicting classifications of pathogenicity", ["Condition A"])])
        ]
        text = assemble_brief(db(), records, [], [], None, clinvar_retrieved=True)
        self.assertIn("Disease relevance: unresolved.", text)
        self.assertIn("Conflicting classifications of pathogenicity", text)

    def test_uniprot_is_function_only_and_only_without_significance(self) -> None:
        hit = UniProtHit("P00000", "Binds a receptor and does not describe a disease.")
        empty = record(classifications=[classification("not provided")], consequences=["missense variant"])
        without = assemble_brief(db(), [empty], [], [], hit, clinvar_retrieved=True)
        self.assertIn("UniProt P00000", without)
        self.assertIn("no clinical significance was found (ClinVar).", without)
        disease = without.split("Disease relevance:")[1]
        self.assertNotIn("UniProt", disease)
        self.assertNotIn("Binds a receptor", disease)
        with_sig = assemble_brief(db(), [record()], [], [], hit, clinvar_retrieved=True)
        self.assertNotIn("UniProt", with_sig)
        self.assertNotIn("Binds a receptor", with_sig)

    def test_no_change_record_is_not_mixed_into_the_alt(self) -> None:
        records = [
            record(accession="VCVDRUG", classifications=[classification("drug response", ["Drug toxicity"])]),
            record(
                accession="VCVREF",
                title="NM_1.2(GENE):c.1=",
                cdna="c.1=",
                consequences=["no sequence alteration"],
                classifications=[classification("Conflicting classifications of pathogenicity")],
                no_change=True,
            ),
        ]
        text = assemble_brief(db(), records, [], [], None, clinvar_retrieved=True)
        self.assertIn("drug response", text)
        self.assertIn("VCVREF", text)
        self.assertNotIn("Disease relevance: unresolved", text)

    def test_brief_does_not_diagnose(self) -> None:
        text = assemble_brief(
            db(),
            [record()],
            [variant()],
            [Paper("Smith J, Doe A", "2020", "A research article", "123")],
            None,
            clinvar_retrieved=True,
        )
        lowered = text.lower()
        self.assertNotIn("diagnos", lowered)
        self.assertNotIn("patient", lowered)
        self.assertNotIn("you have", lowered)
        self.assertIn("PMID:123", text)
        self.assertIn("(Europe PMC).", text)
        self.assertIn("(dbSNP)", text)

    def test_coordinate_conflict_is_kept(self) -> None:
        other = variant()
        other.pos = 999
        text = assemble_brief(db(), [], [other], [], None, clinvar_retrieved=True)
        self.assertIn("Coordinate: unresolved.", text)
        self.assertIn("chr1:100", text)
        self.assertIn("chr1:999", text)


def _summary(accession: str, dbsnp: str, title: str, cdna: str, significance: str, consequences: list[str]) -> dict:
    return {
        "accession": accession,
        "title": title,
        "genes": [{"symbol": "GENE"}],
        "molecular_consequence_list": consequences,
        "germline_classification": {
            "description": significance,
            "review_status": "criteria provided",
            "trait_set": [{"trait_name": "Example"}],
        },
        "variation_set": [
            {
                "cdna_change": cdna,
                "variation_xrefs": [{"db_source": "dbSNP", "db_id": dbsnp}],
            }
        ],
    }


if __name__ == "__main__":
    unittest.main()
