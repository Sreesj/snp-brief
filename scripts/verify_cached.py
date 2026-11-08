"""Independent cross-check of every number the app shows for the cached SNPs.

Run from the repo root: PYTHONPATH=src python3 scripts/verify_cached.py
It recomputes allele frequencies, percentages and carrier rates from the raw cached gnomAD
responses, checks gene and protein change for well-known variants, and flags ACMG frequency
codes that contradict ClinVar without a warning.
"""
import json, re, sys, os, contextlib, io
from pathlib import Path
os.environ["SNP_BRIEF_EXTRAS"] = "on"
from snp_brief.pipeline import gather
from snp_brief.summary import summarize, _plan
from snp_brief.interpret import genotype_rates, acmg_draft
from snp_brief.demo import theme

ROOT = Path("."); CACHE = ROOT / "cache"
ids = [l.strip() for l in open("fixture/snps.txt") if l.strip()]
problems = []
KNOWN = {
    "rs334": ("HBB", "Glu7Val"), "rs6025": ("F5", "Arg534Gln"), "rs429358": ("APOE", "Cys130Arg"),
    "rs7412": ("APOE", "Arg176Cys"), "rs1801133": ("MTHFR", "Ala222Val"), "rs1426654": ("SLC24A5", "Thr111Ala"),
    "rs1799963": ("F2", ""), "rs113993960": ("CFTR", "Phe508del"), "rs28934578": ("TP53", "Arg175His"),
    "rs17822931": ("ABCC11", "Gly180Arg"), "rs3827760": ("EDAR", "Val370Ala"), "rs713598": ("TAS2R38", ""),
    "rs12913832": ("HERC2", ""), "rs4988235": ("MCM6", ""), "rs121913529": ("KRAS", "Gly12"),
}
def bad(rsid, msg): problems.append(f"{rsid}: {msg}")

for rsid in ids:
    if not (CACHE / rsid).is_dir():
        print(f"{rsid:12} not cached here; skipped")
        continue
    with contextlib.redirect_stderr(io.StringIO()):
        item = gather(rsid, cache_dir=CACHE, refresh=False, offline=True)
    raw_dir = CACHE / item.rsid
    # 1. gnomAD counts vs raw
    raw = json.load(open(raw_dir / "gnomad.json")) if (raw_dir / "gnomad.json").exists() else None
    pcts = set()
    raw_counts = {}
    if raw:
        for v in raw["response"].get("variants", []):
            var = v["totals"]["data"]["variant"]
            for cs in ("exome", "genome"):
                d = var.get(cs)
                if d and d.get("an"):
                    raw_counts[(var["variantId"], cs, "all")] = (d["ac"], d["an"])
                    if abs(d["ac"] / d["an"] - (d.get("af") or 0)) > 1e-9:
                        bad(rsid, f"gnomAD af {d.get('af')} != ac/an for {cs}")
            for cs, block in (v.get("populations") or {}).items():
                pv = ((block or {}).get("data") or {}).get("variant") or {}
                for pop in ((pv.get(cs) or {}).get("populations") or []):
                    if pop["an"]:
                        raw_counts[(pv.get("variantId"), cs, pop["id"])] = (pop["ac"], pop["an"])
        for (vid, cs, pid), (ac, an) in raw_counts.items():
            p = ac / an * 100
            for sig in (1, 2, 3):
                pcts.add(float(f"{p:.{sig}g}"))
            pcts.add(round(p, 1)); pcts.add(round(p, 2))
    # parsed callsets equal raw
    for variant in item.variants:
        for cs in variant.callsets:
            key = "exome" if "exome" in cs.name else "genome"
            r = raw_counts.get((variant.variant_id, key, "all"))
            if r and (cs.ac, cs.an) != r:
                bad(rsid, f"parsed {cs.name} {cs.ac}/{cs.an} != raw {r}")
            for pop in cs.populations:
                r = raw_counts.get((variant.variant_id, key, pop.id))
                if r and (pop.ac, pop.an) != r:
                    bad(rsid, f"parsed pop {pop.id} {pop.ac}/{pop.an} != raw {r}")
    # 2. every % in summary and frequency card comes from raw counts
    summary = summarize(item)
    plan = _plan(item)
    texts = summary.text + " " + theme.facts_grid(plan)
    for m in re.finditer(r"(\d+(?:\.\d+)?)%", texts):
        val = float(m.group(1))
        if raw and val not in pcts and not any(abs(val - p) < 1e-9 for p in pcts):
            ctx = texts[max(0, m.start() - 60): m.end() + 10].replace("\n", " ")
            bad(rsid, f"percentage {val}% not derivable from raw gnomAD counts: …{ctx}…")
    # 3. carrier rates recomputed
    for r in genotype_rates(item):
        if abs(r.carriers - (1 - (1 - r.af) ** 2)) > 1e-12 or abs(r.homozygous - r.af ** 2) > 1e-12:
            bad(rsid, f"carrier math wrong for {r.group}")
        if not any(abs(r.af - ac / an) < 1e-12 for ac, an in raw_counts.values()):
            bad(rsid, f"carrier AF {r.af} for {r.group} is not a raw ac/an")
    # 4. ACMG vs ClinVar
    codes = [(c.code, c.strength) for c in acmg_draft(item)]
    head = plan.headline or ""
    benign_freq = any(c in {"BA1", "BS1"} for c, _ in codes)
    if benign_freq and "athogenic" in head and "onflict" not in head and not any(c == "ClinVar" for c, _ in codes):
        bad(rsid, f"ACMG {codes} contradicts ClinVar '{head}' (needs a warning)")
    # 5. dbSNP vs gnomAD position
    if item.db.pos and item.variants and all(abs(v.pos - item.db.pos) > 5 for v in item.variants):
        bad(rsid, f"dbSNP pos {item.db.pos} far from gnomAD {[v.pos for v in item.variants]}")
    # 6. well-known facts (textbook values) for famous variants
    want = KNOWN.get(rsid)
    if want:
        gene, prot = want
        if gene not in (plan.gene or "") + " " + (item.db.genes and ",".join(item.db.genes) or ""):
            bad(rsid, f"gene {plan.gene!r} / dbSNP {item.db.genes} but expected {gene}")
        alltext = summary.text + " " + item.text
        if prot and prot not in alltext:
            bad(rsid, f"expected protein change {prot} not mentioned")
    top = sorted(((ac / an, pid, cs) for (vid, cs, pid), (ac, an) in raw_counts.items() if an >= 2000 and "_" not in pid and pid not in ("XX", "XY")), reverse=True)[:2]
    print(f"{rsid:12} top AF {[(pid, cs, round(af*100, 2)) for af, pid, cs in top]}")
    print(f"{rsid:12} gene={plan.gene!s:12} cons={plan.consequence!s:28} clin={head[:38]!s:38} acmg={codes}")
print()
print("\n".join(problems) or "no problems found")
