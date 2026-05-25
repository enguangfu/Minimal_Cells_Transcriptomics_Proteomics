"""
Annotate the syn3A proteome with Secondary + Tertiary function hierarchy.

STEP 1 (this script): metabolism enzymes.
Iterate the metabolism subsystems in syn3A_metabolism_kinetic_rates.xlsx, collect
every enzyme, resolve it to a syn3A locus tag, and assign a (Secondary, Tertiary)
function from the controlled vocabulary.

Controlled vocabulary: function_hierachy.tsv (forward-filled Primary/Secondary,
one tertiary per row). The classification logic still lives in the dicts/keyword
rules below, but every emitted (Secondary, Tertiary) is validated against this
file — any pair not in it is flagged ILLEGAL_VOCAB in the Review Flag column.

Target proteome table: the 'Comparative Proteomics' sheet of
syn3A_proteome_cplxformation_paper.xlsx (a real .xlsx). It is reshaped to the
canonical Proteome layout — the cross-species / extra localization columns are
dropped and empty Secondary/Tertiary Function columns are inserted — so the
output columns are identical regardless of source workbook.

Subsystems used: Central, Nucleotide, Lipid, Cofactor, Transport, tRNA Charging.
Ignored (no clean per-enzyme loci / per request): SSU Assembly, LSU Assembly,
Gene Expression, Non-Random-Binding Reactions.

Enzyme tokens live in the `Value` column where `Parameter Type == 'Eff Enzyme
Count'` (reaction sheets) or `== 'synthetase'` (tRNA Charging). A token is one of:
  - P_NNNN                -> locus JCVISYN3A_NNNN
  - P_NNNN-P_MMMM / Name-Name  -> hyphen-joined parts, each resolved individually
  - a complex name        -> expanded to member loci via the 'Complexes' sheet
                             ('Genes Products', ';'-separated locus numbers)
  - 'default'             -> skipped

(Secondary, Tertiary) assignment:
  - Nucleotide  -> Central Carbon Metabolism / Nucleotide metabolism
  - Lipid       -> Biosynthesis / Lipid metabolism
  - Cofactor    -> Biosynthesis / Cofactor biosynthesis
  - tRNA Charging -> Translation / tRNA loading
  - Central     -> Central Carbon Metabolism / <tertiary by gene name + product>
                   (Glycolysis / Pentose phosphate metabolism / Pyruvate metabolism
                    / Other central metabolism enzymes; TCA & Carbohydrate excluded
                    as they are absent in syn3A)
  - Transport   -> <by gene name + product>: usually Membrane Transport /
                   {Nucleic acids, Amino acid and peptides, Polyamines, Cofactor,
                    Inorganic, Other Transport}; ATP synthase (atp*) is routed to
                   Energy Metabolism / Oxidative phosphorylation.

Output: syn3A_proteome_tertiary_annotated.tsv — the original Proteome sheet
columns, with Secondary/Tertiary Function filled ONLY for the metabolism enzymes
found here; all other proteins left blank (handled in later steps).
"""

import re
import pandas as pd

META = "Syn3A_annotation/syn3A_metabolism_kinetic_rates.xlsx"
CPLX = "Syn3A_annotation/syn3A_cplx_formation.xlsx"
PROT = "Syn3A_annotation/syn3A_proteome_cplxformation_paper.xlsx"   # 'Comparative Proteomics' sheet
HIER = "Syn3A_annotation/function_hierachy.tsv"                     # controlled (Secondary,Tertiary) vocab
OUT  = "syn3A_proteome_tertiary_annotated.tsv"        # steps 1+2 (curated only)
OUT_AI = "syn3A_proteome_fully_annotated_AI.tsv"      # steps 1+2 + AI best-effort
OUT_REPORT = "annotate_tertiary_function_report.txt"  # run summary

# Subsystems with a fixed (Secondary, Tertiary); None => decide by gene/product.
FIXED = {
    "Nucleotide":    ("Central Carbon Metabolism", "Nucleotide metabolism"),
    "Lipid":         ("Biosynthesis", "Lipid metabolism"),
    "Cofactor":      ("Biosynthesis", "Cofactor biosynthesis"),
    "tRNA Charging": ("Translation", "tRNA loading"),
    "Central":       None,
    "Transport":     None,
}
# Secondary -> implied Primary (to flag Primary/Secondary mismatches for review).
PRIMARY_OF_SECONDARY = {
    "Transcription": "Genetic Information Processing",
    "Translation": "Genetic Information Processing",
    "Folding, Sorting and Degradation": "Genetic Information Processing",
    "DNA Maintenance": "Genetic Information Processing",
    "Signal Transduction": "Environmental Information Processing",
    "Cytoskeleton": "Cellular Processes",
    "Cell Motility": "Cellular Processes",
    "Cell Growth": "Cellular Processes",
    "Defense": "Cellular Processes",
    "Membrane Transport": "Metabolism",
    "Central Carbon Metabolism": "Metabolism",
    "Energy Metabolism": "Metabolism",
    "Biosynthesis": "Metabolism",
    "Other Enzymes": "Metabolism",
    "Drug resistance": "Human Diseases",
    "Infectious Diseases": "Human Diseases",
    "Kegg ortholog defined": "Unclear",
    "No Kegg ortholog": "Unclear",
}
# Enzyme-bearing Parameter Type per sheet.
ENZ_PARAM = {"tRNA Charging": "synthetase"}   # default for the rest: 'Eff Enzyme Count'
# Process order; first subsystem to claim a locus wins (others reported as
# conflicts). Central/Transport (gene/product-resolved) take priority over the
# coarse fixed subsystems, since core carbon enzymes (e.g. pyk, pgk) also appear
# in the Nucleotide/Lipid reaction sheets but belong to Central metabolism.
ORDER = ["tRNA Charging", "Central", "Transport", "Nucleotide", "Lipid", "Cofactor"]


# ── Complex resolution ───────────────────────────────────────────────────────
_cx = pd.read_excel(CPLX, sheet_name="Complexes")
CXMAP = {str(n).strip(): [g.strip().zfill(4) for g in str(gp).split(";") if g.strip()]
         for n, gp in zip(_cx["Name"], _cx["Genes Products"]) if pd.notna(n)}


def resolve(token):
    """Return (loci, unresolved_parts) for one enzyme token."""
    loci, unres = [], []
    tok = str(token).strip()
    if tok.lower() in ("default", "nan", ""):
        return loci, unres
    for part in tok.split("-"):
        part = part.strip()
        m = re.fullmatch(r"P_(\d+)", part)
        if m:
            loci.append(m.group(1).zfill(4))
        elif part in CXMAP:
            loci += CXMAP[part]
        elif re.fullmatch(r"\d+", part):
            loci.append(part.zfill(4))
        else:
            unres.append(part)
    return loci, unres


# ── Gene/product-based tertiary classifiers (Central & Transport) ────────────
def _has(text, kws):
    t = text.lower()
    return any(k in t for k in kws)


def central_tertiary(gn, gp):
    s = f"{gn} {gp}"
    if _has(s, ["phosphofructokinase", "bisphosphate aldolase", "triose-phosphate isomerase",
                "triosephosphate isomerase", "glyceraldehyde-3-phosphate dehydrogenase",
                "phosphoglycerate kinase", "phosphoglycerate mutase", "phosphopyruvate hydratase",
                "enolase", "pyruvate kinase", "glucose-6-phosphate isomerase"]):
        return "Glycolysis"
    if _has(s, ["ribulose-phosphate", "transketolase", "transaldolase", "ribose 5-phosphate isomerase",
                "ribose-5-phosphate isomerase", "deoxyribose-phosphate aldolase", "phosphopentomutase",
                "phosphoribosylpyrophosphate"]):
        return "Pentose phosphate metabolism"
    if _has(s, ["lactate dehydrogenase", "pyruvate dehydrogenase", "dihydrolipoyl dehydrogenase",
                "acetate kinase", "acetyltransferase", "alpha-keto acid dehydrogenase"]):
        return "Pyruvate metabolism"
    return "Other central metabolism enzymes"


def transport_class(gn, gp):
    s = f"{gn} {gp}"
    if re.match(r"atp[a-z]", str(gn).strip().lower()) or _has(s, ["atp synthase", "f1 atp", "f0 atp"]):
        return ("Energy Metabolism", "Oxidative phosphorylation")
    if _has(s, ["nucleoside", "nucleotide", "purine", "pyrimidine"]):
        return ("Membrane Transport", "Nucleic acids")
    if _has(s, ["spermidine", "putrescine", "polyamine"]):
        return ("Membrane Transport", "Polyamines")
    if _has(s, ["amino acid", "oligopeptide", "peptide", "glutamate", "glutamine"]):
        return ("Membrane Transport", "Amino acid and peptides")
    if _has(s, ["thiamine", "riboflavin", "folate", "flavin", "cobalamin", "biotin",
                "pantothenate", "coa", "ecf transporter", "ecf "]):
        return ("Membrane Transport", "Cofactor")
    if _has(s, ["phosphate", "potassium", "sodium", "magnesium", "metal", "zinc",
                "manganese", "iron", "cation", "chloride"]):
        return ("Membrane Transport", "Inorganic")
    return ("Membrane Transport", "Other Transport")


# ── AI knowledge classifier for the remaining (unmodeled) proteins ───────────
# Best-effort (Secondary, Tertiary) from gene name + product. Returns (None,None)
# when no rule matches. Ordered specific -> generic. For REVIEW by the user.
_ACTIVITY = ["transferase", "hydrolase", "kinase", "reductase", "oxidase", "synthase",
             "synthetase", "methyltransferase", "phosphohydrolase", "nuclease",
             "peptidase", "protease", "atpase", "isomerase", "epimerase", "ligase",
             "phosphatase", "dehydrogenase", "deaminase", "lipoprotein", "permease"]


def ai_classify(gn, gp, primary):
    s = f"{gn} {gp}".lower()

    def h(*kw):
        return any(k in s for k in kw)

    # Translation: tRNA aminoacylation / modification / processing
    if h("trna ligase", "aminoacyl-trna", "trna synthetase", "--trna ligase"):
        return ("Translation", "tRNA loading")
    if h("ssra-binding"):
        return ("Translation", "Translation factors")
    if "trna" in s and h("methyltransferase", "pseudouridine", "synthase", "thiouridine",
                         "dihydrouridine", "threonylcarbamoyl", "lysidine",
                         "carboxymethylaminomethyl", "amidotransferase", "4-thiouridine"):
        return ("Translation", "tRNA biogenesis")
    # Ribosome biogenesis (rRNA processing / modification / assembly factors)
    if h("16s", "23s", "5s rrna", "rrna", "ribosome biogenesis", "ribosome assembly",
         "ribosome-binding factor", "ribosome binding factor", "ribosome small subunit",
         "ribosome gtpase", "pre-16s", "ribonuclease m5", "ribonuclease iii",
         "rrna maturation", "l16-binding", "ribosomal protein l27", "rimm", "rimp",
         "50s subunit-maturation", "small subunit-dependent gtpase"):
        return ("Translation", "Ribosome biogenesis")
    if h("ribosomal protein") and not h("maturation"):
        return ("Translation", "Ribosome")
    # Translation factors
    if h("elongation factor", "initiation factor", "release factor", "ribosome recycling",
         "peptide chain release", "peptide deformylase", "peptidyl-trna hydrolase",
         "aminoacyl-trna hydrolase", "elongation factor p"):
        return ("Translation", "Translation factors")
    # Transcription
    if h("rna polymerase", "sigma factor"):
        return ("Transcription", "Transcription machinery")
    if h("antitermination", "transcriptional regulator", "transcriptional repressor",
         "transcription factor", "transcription elongation", "transcription termination",
         "transcription antitermination"):
        return ("Transcription", "Transcription factors")
    # Folding, Sorting and Degradation
    if h("chaperone", "chaperonin", "heat shock", "nucleotide exchange factor",
         "trigger factor", "clp protease subunit b"):
        return ("Folding, Sorting and Degradation", "Chaperones and folding catalysts")
    if h("translocase", "signal recognition particle", "membrane protein insertase",
         "signal peptidase", "preprotein", "srp"):
        return ("Folding, Sorting and Degradation", "Protein export")
    if h("ribonuclease r", "ribonuclease p", "degradosome", "exoribonuclease"):
        return ("Folding, Sorting and Degradation", "RNA degradation")
    if h("peptidase", "protease", "endopeptidase", "oligopeptidase", "prolidase",
         "aminopeptidase", "deglycase"):
        return ("Folding, Sorting and Degradation", "Peptidases")
    # DNA Maintenance
    if h("dna polymerase iii"):
        return ("DNA Maintenance", "DNA replication complex")
    if h("replication initiator", "replication initiation", "primosom"):
        return ("DNA Maintenance", "DNA replication control")
    if h("dna polymerase", "dna primase", "dna ligase", "dna helicase",
         "single-stranded dna-binding", "flap endonuclease", "replicative dna helicase",
         "rnase hi", "ribonuclease hi"):
        return ("DNA Maintenance", "DNA replication")
    if h("excinuclease", "excision repair"):
        return ("DNA Maintenance", "Nucleotide excision repair")
    if h("glycosylase", "endonuclease iv", "formamidopyrimidine", "deoxyribonuclease iv"):
        return ("DNA Maintenance", "Base excision repair")
    if h("mismatch"):
        return ("DNA Maintenance", "Mismatch repair")
    if h("recombination", "recombinase"):
        return ("DNA Maintenance", "Homologous recombination")
    if h("dna repair"):
        return ("DNA Maintenance", "DNA repair and recombination proteins")
    if h("topoisomerase", "dna-binding protein", "nucleoid", "histone-like", "whia"):
        return ("DNA Maintenance", "Chromosome-related")
    # Cellular Processes
    if h("cell division", "ftsz", "ftsa", "sepf", "gpsb", "divic", "ezra"):
        return ("Cell Growth", "Cell division")
    if h("toxin", "antitoxin", "restriction", "crispr", "defense"):
        return ("Defense", "Prokaryotic defense system")
    # Human Diseases
    if h("tetracycline resistance", "antibiotic resistance", "antimicrobial resistance",
         "ribosomal protection"):
        return ("Drug resistance", "Antimicrobial resistance genes")
    # Metabolism — transport
    if h("efflux"):
        return ("Membrane Transport", "Drug resistance")
    if h("pts ", "phosphotransferase system", "pts sugar", "pts glucose", "phosphocarrier"):
        return ("Membrane Transport", "Carbohydrates")
    if h("abc transporter", "mfs transporter", "permease", "symporter", "transporter",
         "c4-dicarboxylate"):
        if h("phosphate"):
            return ("Membrane Transport", "Inorganic")
        if h("amino acid", "peptide"):
            return ("Membrane Transport", "Amino acid and peptides")
        if h("nucleoside", "nucleotide"):
            return ("Membrane Transport", "Nucleic acids")
        return ("Membrane Transport", "Other Transport")
    if h("phosphate transport"):
        return ("Membrane Transport", "Inorganic")
    # Metabolism — biosynthesis / central
    if h("acyl carrier", "glycolipid synthase", "diacylglyceryl", "fatty acid", "flippase"):
        return ("Biosynthesis", "Lipid metabolism")
    if h("iron-sulfur", "cysteine desulfurase", "adenylyl-sulfate kinase", "sulfate"):
        return ("Biosynthesis", "Cofactor biosynthesis")
    if h("histidine", "methionine adenosyltransferase", "amino acid biosynthesis"):
        return ("Biosynthesis", "Amino acid metabolism")
    if h("nucleosidase", "ntp pyrophosphatase", "pyrophosphohydrolase", "methylthioadenosine"):
        return ("Central Carbon Metabolism", "Nucleotide metabolism")
    if h("serine/threonine protein kinase", "protein kinase"):
        return ("Other Enzymes", "Protein kinases")
    if primary == "Metabolism" and h("mannose", "glucosamine", "acetylmannosamine",
                                     "epimerase", "isomerase", "deaminase"):
        return ("Central Carbon Metabolism", "Other central metabolism enzymes")
    if h("reductase", "thioredoxin", "peroxiredoxin", "oxidase", "dehydrogenase",
         "hydrolase", "transferase", "phosphohydrolase", "phosphatase"):
        return ("Other Enzymes", "Other enzymes")
    return (None, None)


def classify_remaining(primary, gn, gp):
    sec, ter = ai_classify(gn, gp, primary)
    if ter is not None:
        return sec, ter
    s = f"{gn} {gp}".lower()
    if "uncharacterized protein" in s and not any(a in s for a in _ACTIVITY):
        return ("No Kegg ortholog", "Function unknown")
    if any(a in s for a in _ACTIVITY) or "uncharacterized" in s:
        return ("Kegg ortholog defined", "General function prediction only")
    return ("No Kegg ortholog", "Function unknown")


# ── Proteome (target table) ──────────────────────────────────────────────────
# Source: 'Comparative Proteomics' sheet of the paper workbook (a real .xlsx).
# Reshape it to the canonical Proteome layout — drop the cross-species / extra
# localization columns, (re)insert empty Secondary/Tertiary Function — so the
# output columns are identical to the previous tertiaryfunction.xlsx-driven run.
PROTEOME_COLS = ["Locus Tag", "Gene Name", "Gene Product", "Protein Length",
                 "Exp. Ptn Cnt", "Essentiality", "Primary Function",
                 "Secondary Function", "Tertiary Function", "Localization",
                 "Sim. Initial Ptn Cnt"]
prot = pd.read_excel(PROT, sheet_name="Comparative Proteomics")
prot = prot[prot["Locus Tag"].astype(str).str.match(r"JCVISYN3A_\d+", na=False)].copy()
for c in ("Secondary Function", "Tertiary Function"):
    if c not in prot.columns:
        prot[c] = pd.NA
prot = prot[PROTEOME_COLS].copy()
prot["_num"] = prot["Locus Tag"].astype(str).str.extract(r"(\d+)$")[0]
info = prot.set_index("_num")[["Gene Name", "Gene Product"]].to_dict("index")
ORIG_COLS = [c for c in prot.columns if c != "_num"]


# ── Collect enzymes per subsystem ────────────────────────────────────────────
meta = pd.ExcelFile(META)
assigned = {}          # num -> (subsystem, secondary, tertiary)
conflicts, unresolved, unmatched = [], [], []
per_subsystem = {}

for sub in ORDER:
    df = pd.read_excel(meta, sheet_name=sub)
    ptype = ENZ_PARAM.get(sub, "Eff Enzyme Count")
    tokens = df.loc[df["Parameter Type"].astype(str) == ptype, "Value"].dropna().astype(str)
    loci = set()
    for tok in tokens:
        ll, uu = resolve(tok)
        loci.update(ll)
        for u in uu:
            unresolved.append((sub, tok, u))
    per_subsystem[sub] = len(loci)

    for num in sorted(loci):
        if num not in info:
            unmatched.append((sub, num))
            continue
        if num in assigned:
            conflicts.append((num, assigned[num][0], sub))
            continue
        gn = str(info[num].get("Gene Name") or "")
        gp = str(info[num].get("Gene Product") or "")
        if FIXED[sub] is not None:
            sec, ter = FIXED[sub]
        elif sub == "Central":
            sec, ter = "Central Carbon Metabolism", central_tertiary(gn, gp)
        else:  # Transport
            sec, ter = transport_class(gn, gp)
        assigned[num] = (sub, sec, ter)


# ── STEP 2: complex-based annotation (all Genetic Information Processing) ─────
# Expand the named complexes to member loci (Complexes 'Genes Products') and
# assign (Secondary, Tertiary). ExoVII left out (function unknown per user).
COMPLEX_FUNC = {
    "RNAP":        ("Transcription", "Transcription machinery"),
    "TopoIV":      ("DNA Maintenance", "Chromosome-related"),
    "Gyrase":      ("DNA Maintenance", "Chromosome-related"),
    "SMC":         ("DNA Maintenance", "Chromosome-related"),
    "ExoVII":      ("DNA Maintenance", "DNA repair and recombination proteins"),
    "SecYEGDF":    ("Folding, Sorting and Degradation", "Protein export"),
    "Ribosome":    ("Translation", "Ribosome"),
    "Degradosome": ("Folding, Sorting and Degradation", "RNA degradation"),
}
cplx_counts = {}
for cxname, (sec, ter) in COMPLEX_FUNC.items():
    n_set = 0
    for num in CXMAP.get(cxname, []):
        if num not in info:
            unmatched.append((f"cplx:{cxname}", num))
            continue
        if num in assigned:
            conflicts.append((num, assigned[num][0], f"cplx:{cxname}"))
            continue
        assigned[num] = (f"cplx:{cxname}", sec, ter)
        n_set += 1
    cplx_counts[cxname] = n_set


# ── Controlled vocabulary (single source of truth for legal functions) ───────
# function_hierachy.tsv is forward-filled (Primary/Secondary repeat down as
# blanks); one tertiary per row after the cleanup. Any assigned (Secondary,
# Tertiary) not present here is flagged ILLEGAL_VOCAB for review.
_hier = pd.read_csv(HIER, sep="\t")
_hier["Secondary"] = _hier["Secondary"].ffill()
LEGAL_PAIRS = {(str(s).strip(), str(t).strip())
               for s, t in zip(_hier["Secondary"], _hier["Tertiary"])
               if pd.notna(t) and str(t).strip()}


# ── Review Flag helper (ambiguity for manual inspection) ─────────────────────
conflict_other = {num: also for num, first, also in conflicts}


def review_flag(num, secondary, tertiary, sheet_primary, is_ai):
    flags = []
    if num in conflict_other:
        flags.append(f"CONFLICT(also:{conflict_other[num]})")
    if pd.notna(secondary):
        implied = PRIMARY_OF_SECONDARY.get(secondary)
        if implied and str(sheet_primary) not in ("", "nan", "None") and implied != str(sheet_primary):
            flags.append("PRIMARY_MISMATCH")
    if pd.notna(secondary) and pd.notna(tertiary) and \
            (str(secondary).strip(), str(tertiary).strip()) not in LEGAL_PAIRS:
        flags.append("ILLEGAL_VOCAB")
    if is_ai:
        flags.append("AI")
    return "; ".join(flags)


OUT_COLS = ORIG_COLS + ["Review Flag"]

# ── Fill the proteome copy and write TSV (steps 1+2 only) ────────────────────
prot["Secondary Function"] = prot["_num"].map(lambda n: assigned[n][1] if n in assigned else pd.NA)
prot["Tertiary Function"]  = prot["_num"].map(lambda n: assigned[n][2] if n in assigned else pd.NA)
prot["Review Flag"] = [review_flag(r["_num"], r["Secondary Function"], r["Tertiary Function"],
                                   r["Primary Function"], False)
                       for _, r in prot.iterrows()]
prot[OUT_COLS].to_csv(OUT, sep="\t", index=False)


# ── STEP 3: AI knowledge classification of all remaining proteins ────────────
# Best-effort; keeps the existing Primary Function, fills Secondary/Tertiary.
ai = prot.copy()
n_ai = 0
for idx, r in ai.iterrows():
    if pd.notna(r["Tertiary Function"]):
        continue
    sec, ter = classify_remaining(str(r.get("Primary Function") or ""),
                                  str(r.get("Gene Name") or ""), str(r.get("Gene Product") or ""))
    ai.at[idx, "Secondary Function"] = sec
    ai.at[idx, "Tertiary Function"] = ter
    n_ai += 1
ai["Review Flag"] = [review_flag(r["_num"], r["Secondary Function"], r["Tertiary Function"],
                                 r["Primary Function"], r["_num"] not in assigned)
                     for _, r in ai.iterrows()]
ai[OUT_COLS].to_csv(OUT_AI, sep="\t", index=False)
n_blank_ai = int(ai["Tertiary Function"].isna().sum())


# ── Report (printed to console AND written to OUT_REPORT) ────────────────────
_R = []


def say(line=""):
    print(line)
    _R.append(line)


n_metab = (prot["Primary Function"].astype(str) == "Metabolism").sum()
say("# syn3A proteome tertiary-function annotation report")
say(f"Proteome rows (valid loci): {len(prot)}; Primary=='Metabolism': {n_metab}")
say("\nStep 1 — enzyme loci per metabolism subsystem:")
for sub in ORDER:
    say(f"  {sub:<14} {per_subsystem[sub]}")
say("Step 2 — complex members assigned:")
for cx, n in cplx_counts.items():
    say(f"  {cx:<14} {n}")
say(f"\nCurated total (steps 1+2): {len(assigned)}")
say(f"Conflicts (locus claimed by >1 subsystem, kept first): {len(conflicts)}")
for num, first, also in conflicts:
    say(f"  {num}: {first} (kept) vs {also}")
say(f"Unresolved tokens: {len(unresolved)}")
for sub, tok, u in unresolved:
    say(f"  [{sub}] {tok} -> {u}")
say(f"Loci not in proteome: {len(unmatched)}")
for sub, num in unmatched:
    say(f"  [{sub}] {num}")

say("\n--- Central / Transport tertiary (gene-based) for review ---")
rev = prot[prot["_num"].isin([n for n, v in assigned.items() if v[0] in ("Central", "Transport")])]
for _, r in rev.sort_values("_num").iterrows():
    say(f"  {r['Locus Tag']}  {str(r['Gene Name']):<10} {str(r['Secondary Function']):<26} "
        f"{str(r['Tertiary Function']):<28} | {str(r['Gene Product'])[:45]}")

n_unannot = int(prot["Tertiary Function"].isna().sum())
n_annot = len(prot) - n_unannot
say(f"\nCurated (steps 1+2) tertiary annotated: {n_annot} / {len(prot)}  |  unannotated: {n_unannot}")
say(f"STEP 3 AI-classified: {n_ai} proteins  |  still blank after AI: {n_blank_ai}")
_aifilled = ai[~ai["_num"].isin(list(assigned))]
say(f"\nStep 3 AI tertiary distribution ({n_ai} AI-filled rows):")
say(_aifilled["Tertiary Function"].value_counts(dropna=False).to_string())

say("\nReview Flag counts (in fully-annotated AI table):")
nz = ai["Review Flag"].astype(str).str.len() > 0
say(f"  CONFLICT:         {int(ai['Review Flag'].str.contains('CONFLICT').sum())}")
say(f"  PRIMARY_MISMATCH: {int(ai['Review Flag'].str.contains('PRIMARY_MISMATCH').sum())}")
say(f"  ILLEGAL_VOCAB:    {int(ai['Review Flag'].str.contains('ILLEGAL_VOCAB').sum())}")
say(f"  AI:               {int(ai['Review Flag'].str.contains('AI').sum())}")
say(f"Controlled vocab: {len(LEGAL_PAIRS)} legal (Secondary,Tertiary) pairs from {HIER}")
_illegal = ai[ai["Review Flag"].str.contains("ILLEGAL_VOCAB")]
if len(_illegal):
    say("  ILLEGAL_VOCAB pairs (not in function_hierachy.tsv):")
    for _, r in _illegal.iterrows():
        say(f"    {r['Locus Tag']}  {r['Secondary Function']} / {r['Tertiary Function']}")
say("  PRIMARY_MISMATCH loci (Primary vs implied):")
for _, r in ai[ai["Review Flag"].str.contains("PRIMARY_MISMATCH")].iterrows():
    say(f"    {r['Locus Tag']}  {str(r['Gene Name']):<10} sheet={r['Primary Function']:<30} "
        f"-> {r['Secondary Function']} / {r['Tertiary Function']}")

say("\n--- Step 3 AI calls, by Primary Function (for one-by-one review) ---")
for prim in ["Genetic Information Processing", "Metabolism", "Cellular Processes",
             "Human Diseases", "Unclear"]:
    block = _aifilled[_aifilled["Primary Function"].astype(str) == prim]
    if block.empty:
        continue
    say(f"\n### {prim}  (n={len(block)})")
    for _, r in block.sort_values("Locus Tag").iterrows():
        say(f"  {r['Locus Tag']}  {str(r['Gene Name']):<12} {str(r['Secondary Function']):<34} "
            f"{str(r['Tertiary Function']):<34} | {str(r['Gene Product'])[:45]}")

say(f"\nSaved: {OUT}")
say(f"Saved: {OUT_AI}")
say(f"Saved: {OUT_REPORT}")
with open(OUT_REPORT, "w") as fh:
    fh.write("\n".join(_R) + "\n")
