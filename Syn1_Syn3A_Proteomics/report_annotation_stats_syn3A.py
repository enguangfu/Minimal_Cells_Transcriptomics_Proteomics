"""
Report statistics for the curated syn3A proteome functional annotation.

INPUT
  syn3A_proteome_fully_annotated_Revised.xlsx (sheet 'Proteome') — 455 proteins
  with a hand-curated Primary / Secondary / Tertiary function hierarchy. Produced
  by annotate_tertiary_function_syn3A.py (steps 1-3) and then manually reviewed
  (CONFLICT / PRIMARY_MISMATCH / AI rows corrected; `Mismatch Solved` marks the
  PRIMARY_MISMATCH rows that were adjudicated). The .xlsx is the curation master
  (keeps color highlights / sort state that a TSV cannot). Controlled vocabulary:
  Syn3A_annotation/function_hierachy.tsv.

WHAT IT REPORTS  (printed to console AND written to OUT_REPORT)
  1. Overview          — N proteins, essentiality, localization.
  2. Functional tree   — Primary > Secondary > Tertiary counts (nested).
  3. Abundance shares  — fraction of the *measured* proteome (Exp. Ptn Cnt copies)
                         carried by each Primary / Secondary category, next to its
                         gene-count share (a few classes dominate copy number even
                         with modest gene counts — e.g. ribosomal proteins).
  4. Cross-tabs        — Primary x Essentiality, Primary x Localization.
  5. Vocab validation  — any (Secondary, Tertiary) pair absent from
                         function_hierachy.tsv (flags manual-edit slips).
  6. Curation trail    — Review Flag breakdown + Mismatch Solved count.

Protein sequences are attached from FASTA (Genomes_Input/syn3A_ptns.fasta, keyed
by locus tag; trailing '*' stop removed) and exported in both the xlsx and html.

ALSO WRITES
  OUT_TABLE — the full Primary>Secondary>Tertiary table with n_genes, summed
  Exp. Ptn Cnt, and both shares (% of proteins, % of copies); machine-readable.
  OUT_XLSX  — derived annotated workbook: columns reordered per spec with Protein
  Sequence after Gene Product (curation columns appended last). The _Revised
  curation master is the read-only INPUT and is never modified.
  OUT_FIG   — grouped horizontal-bar figure of every tertiary function, bars
  grouped under a bold Primary header and colored by Primary function.
  OUT_HTML  — self-contained interactive version: clickable composition bars;
  clicking a tertiary function filters the protein table below (sortable headers,
  CSV/Excel download of the current view, width-clamped expandable sequence cell).
"""

import json
from html import escape as esc

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Match the HTML's system sans-serif look (falls back to DejaVu Sans on Linux).
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"]})

IN         = "syn3A_proteome_fully_annotated_Revised.xlsx"
IN_SHEET   = 0     # first sheet (robust to sheet renames in the curation master)
FASTA      = "../Genomes_Input/syn3A_ptns.fasta"   # protein sequences (locus-tag headers)
PROT2026   = "syn3a_proteomics_summary_2026.csv"   # has copy_number_2026 per locus_tag
RNA_ABUND  = "../Syn3A_Corr_RNA_Proteins/syn3A_rna_abundances.tsv"  # mRNA copies/cell + Illumina TPM per locus (Calc_Abundances.py)
HIER       = "Syn3A_annotation/function_hierachy.tsv"
OUT_REPORT = "syn3A_annotation_stats.txt"
OUT_TABLE  = "syn3A_function_hierarchy_counts.tsv"
OUT_FIG    = "syn3A_tertiary_function_composition.pdf"
OUT_HTML   = "syn3A_tertiary_function_composition.html"
OUT_XLSX   = "syn3A_proteome_annotated.xlsx"        # derived: reordered + Protein Sequence

# Final column order for the annotated xlsx + HTML table (abundance counts after
# Gene Product; Protein Sequence/Length after the function levels; curation
# columns appended last). Missing columns are skipped. The loaded "Exp. Ptn Cnt"
# is the 2019 measurement -> renamed "Exp. Ptn Cnt 2019"; "Exp. Ptn Cnt 2026"
# (copy_number_2026 from PROT2026, rounded up) is inserted right after it.
TABLE_COLS = ["Locus Tag", "Gene Name", "Gene Product", "Exp. Ptn Cnt 2019",
              "Exp. Ptn Cnt 2026", "Sim. Initial Ptn Cnt",
              "mRNA Copies/Cell", "Illumina TPM", "Localization",
              "Essentiality", "Primary Function", "Secondary Function",
              "Tertiary Function", "Protein Sequence", "Protein Length",
              "Review Flag", "Mismatch Solved"]
SEQCOL = "Protein Sequence"
ABUND      = "Exp. Ptn Cnt 2019"   # measured protein copy number per cell (2019)
PRIM, SEC, TER = "Primary Function", "Secondary Function", "Tertiary Function"

# Stable Primary-Function palette (muted, print-friendly).
PRIM_COLORS = {
    "Genetic Information Processing":       "#3b6db3",  # blue
    "Metabolism":                          "#3f9e5a",  # green
    "Unclear":                             "#9aa0a6",  # grey
    "Cellular Processes":                  "#8e6bb1",  # purple
    "Environmental Information Processing": "#2aa6a0",  # teal
    "Exogenous":                           "#c0654e",  # terracotta
}


# ── Load ─────────────────────────────────────────────────────────────────────
def read_fasta(path):
    """{header_first_token: sequence} from a FASTA file."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].split()[0], []
            elif ln:
                buf.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


d = pd.read_excel(IN, sheet_name=IN_SHEET)
d = d.rename(columns={"Exp. Ptn Cnt": "Exp. Ptn Cnt 2019"})   # the loaded count is 2019
d[ABUND] = pd.to_numeric(d[ABUND], errors="coerce")
N = len(d)
total_copies = d[ABUND].sum()

# attach protein sequences (strip the trailing '*' stop translation)
_seq = read_fasta(FASTA)
d[SEQCOL] = d["Locus Tag"].map(lambda t: _seq.get(t, "").rstrip("*"))
_missing_seq = int((d[SEQCOL] == "").sum())

# attach 2026 copy numbers (rounded UP to integer; missing -> blank) right after 2019
_cn26 = pd.read_csv(PROT2026).set_index("locus_tag")["copy_number_2026"]
d["Exp. Ptn Cnt 2026"] = np.ceil(d["Locus Tag"].map(_cn26)).astype("Int64")
_missing_2026 = int(d["Exp. Ptn Cnt 2026"].isna().sum())

# attach Syn3A transcription (absolute mRNA copies/cell + Illumina sense TPM) by locus tag
_rna = pd.read_csv(RNA_ABUND, sep="\t").set_index("locus_tag")
d["mRNA Copies/Cell"] = d["Locus Tag"].map(_rna["copies_per_cell"]).round(2)
d["Illumina TPM"] = d["Locus Tag"].map(_rna["Illumina_sense_TPM"]).round(1)
_missing_rna = int(d["mRNA Copies/Cell"].isna().sum())

# resolve the requested column order to those actually present
ORDER_COLS = [c for c in TABLE_COLS if c in d.columns]

# ── Derived annotated workbook (reordered + Protein Sequence) ────────────────
d[ORDER_COLS].to_excel(OUT_XLSX, sheet_name="Proteome", index=False)


# ── Report sink ──────────────────────────────────────────────────────────────
_R = []


def say(line=""):
    print(line)
    _R.append(line)


def share_table(col, title):
    """Per-category gene-count share and copy-number (Exp. Ptn Cnt) share."""
    g = d.groupby(col, dropna=False).agg(n_genes=("Locus Tag", "size"),
                                         copies=(ABUND, "sum"))
    g["pct_genes"] = 100 * g["n_genes"] / N
    g["pct_copies"] = 100 * g["copies"] / total_copies
    g = g.sort_values("copies", ascending=False)
    say(f"\n{title}")
    say(f"  {'category':<38} {'genes':>6} {'%gene':>7} {'copies':>10} {'%copy':>7}")
    for cat, r in g.iterrows():
        say(f"  {str(cat):<38} {int(r['n_genes']):>6} {r['pct_genes']:>6.1f}% "
            f"{int(r['copies']):>10} {r['pct_copies']:>6.1f}%")
    return g


# ── 1. Overview ──────────────────────────────────────────────────────────────
say("# syn3A proteome — curated functional-annotation statistics")
say(f"Input: {IN}")
say(f"Proteins: {N} | total measured copies (Exp. Ptn Cnt 2019): {int(total_copies):,} "
    f"| genes with a name: {d['Gene Name'].notna().sum()}")
say(f"Protein sequences attached from {FASTA}: {N - _missing_seq}/{N}"
    + (f"  (MISSING: {_missing_seq})" if _missing_seq else ""))
say(f"Exp. Ptn Cnt 2026 attached from {PROT2026} (copy_number_2026, rounded up): "
    f"{N - _missing_2026}/{N}" + (f"  (MISSING: {_missing_2026})" if _missing_2026 else ""))

say("\nEssentiality:")
for k, v in d["Essentiality"].value_counts(dropna=False).items():
    say(f"  {str(k):<22} {v:>4}  ({100*v/N:.1f}%)")
say("\nLocalization:")
for k, v in d["Localization"].value_counts(dropna=False).items():
    say(f"  {str(k):<22} {v:>4}  ({100*v/N:.1f}%)")


# ── 2. Functional hierarchy tree ─────────────────────────────────────────────
say("\n" + "=" * 70)
say("Functional hierarchy  (Primary > Secondary > Tertiary; n genes)")
say("=" * 70)
grp = (d.groupby([PRIM, SEC, TER], dropna=False)
         .agg(n_genes=("Locus Tag", "size"), copies=(ABUND, "sum"))
         .reset_index())
for prim in grp[PRIM].drop_duplicates():
    pblk = grp[grp[PRIM] == prim]
    say(f"\n{prim}  (n={int(pblk['n_genes'].sum())})")
    for sec in pblk[SEC].drop_duplicates():
        sblk = pblk[pblk[SEC] == sec]
        say(f"  {sec}  (n={int(sblk['n_genes'].sum())})")
        for _, r in sblk.sort_values("n_genes", ascending=False).iterrows():
            say(f"      {str(r[TER]):<44} {int(r['n_genes']):>4}")


# ── 3. Abundance-weighted shares ─────────────────────────────────────────────
say("\n" + "=" * 70)
say("Composition by gene count vs measured copy number")
say("=" * 70)
share_table(PRIM, "By Primary Function:")
share_table(SEC, "By Secondary Function:")


# ── 4. Cross-tabs ────────────────────────────────────────────────────────────
say("\n" + "=" * 70)
say("Cross-tabulations")
say("=" * 70)
say("\nPrimary Function x Essentiality:")
say(pd.crosstab(d[PRIM], d["Essentiality"], margins=True).to_string())
say("\nPrimary Function x Localization:")
say(pd.crosstab(d[PRIM], d["Localization"], margins=True).to_string())


# ── 5. Vocabulary validation ─────────────────────────────────────────────────
say("\n" + "=" * 70)
say("Controlled-vocabulary validation")
say("=" * 70)
h = pd.read_csv(HIER, sep="\t")
h["Secondary"] = h["Secondary"].ffill()
legal = {(str(s).strip(), str(t).strip())
         for s, t in zip(h["Secondary"], h["Tertiary"])
         if pd.notna(t) and str(t).strip()}
say(f"Legal (Secondary, Tertiary) pairs in {HIER}: {len(legal)}")
illegal = d[[not ((str(s).strip(), str(t).strip()) in legal)
             for s, t in zip(d[SEC], d[TER])]]
if illegal.empty:
    say("All emitted (Secondary, Tertiary) pairs are legal. ✓")
else:
    say(f"ILLEGAL pairs (not in vocab): {len(illegal)} protein(s) — review:")
    for _, r in illegal.iterrows():
        say(f"  {r['Locus Tag']}  {str(r['Gene Name']):<10} "
            f"{r[SEC]} / {r[TER]}")


# ── 6. Curation trail ────────────────────────────────────────────────────────
say("\n" + "=" * 70)
say("Curation provenance")
say("=" * 70)
say("Review Flag breakdown:")
for k, v in d["Review Flag"].value_counts(dropna=False).items():
    say(f"  {str(k):<28} {v:>4}")
if "Mismatch Solved" in d.columns:
    say(f"Mismatch Solved (adjudicated PRIMARY_MISMATCH rows): "
        f"{int(d['Mismatch Solved'].notna().sum())}")


# ── Write machine-readable hierarchy table ───────────────────────────────────
tbl = grp.copy()
tbl["pct_proteins"] = (100 * tbl["n_genes"] / N).round(2)
tbl["pct_copies"] = (100 * tbl["copies"] / total_copies).round(2)
tbl = tbl.rename(columns={"copies": "exp_ptn_copies"})
# tbl.sort_values([PRIM, SEC, "n_genes"], ascending=[True, True, False]).to_csv(
#     OUT_TABLE, sep="\t", index=False)

# ── Figure: tertiary-function composition (3x2 panels, one per Primary) ──────
# In GIP and Metabolism the tertiary bars are clustered under their Secondary
# function (secondaries ordered by total protein count); the four smaller
# Primaries show tertiary bars directly.
SECGROUP = {"Genetic Information Processing", "Metabolism"}


def _panel_rows(prim):
    """Ordered row list for one Primary panel, top-to-bottom.
    rows: ('sec', secondary, total) header  |  ('bar', tertiary, count)."""
    sub = d[d[PRIM] == prim]
    rows = []
    if prim in SECGROUP:                       # cluster tertiary under secondary
        for sec in sub.groupby(SEC).size().sort_values(ascending=False).index:
            ssub = sub[sub[SEC] == sec].groupby(TER).size().sort_values(ascending=False)
            rows.append(("sec", sec, int(ssub.sum())))
            rows += [("bar", t, int(c)) for t, c in ssub.items()]
    else:                                       # tertiary bars only
        rows = [("bar", t, int(c))
                for t, c in sub.groupby(TER).size().sort_values(ascending=False).items()]
    return rows, len(sub)


def _draw_panel(ax, prim, ymax):
    rows, n_prot = _panel_rows(prim)
    color = PRIM_COLORS.get(prim, "#777777")
    xmax = max([r[2] for r in rows if r[0] == "bar"], default=1)

    # alternating faint background band per secondary block
    if prim in SECGROUP:
        seg = [i for i, r in enumerate(rows) if r[0] == "sec"] + [len(rows)]
        for k, (a, b) in enumerate(zip(seg[:-1], seg[1:])):
            ys = [ymax - 1 - i for i in range(a, b)]
            ax.axhspan(min(ys) - 0.5, max(ys) + 0.5, color=color,
                       alpha=0.12 if k % 2 else 0.06, lw=0)

    yticks, ylabels, header_idx = [], [], []
    for i, r in enumerate(rows):
        y = ymax - 1 - i
        yticks.append(y)
        if r[0] == "sec":
            ylabels.append(f"{r[1]}  ({r[2]}, {round(100 * r[2] / N)}%)")
            header_idx.append(i)
        else:
            ylabels.append(r[1])
            ax.barh(y, r[2], height=0.74, color=color, edgecolor="white", linewidth=0.5)
            ax.text(r[2] + xmax * 0.02, y, str(r[2]), va="center", fontsize=7, color="#333")

    ax.set_yticks(yticks)
    lbls = ax.set_yticklabels(ylabels, fontsize=8)
    for k in header_idx:                       # bold + colored secondary headers
        lbls[k].set_fontweight("bold")
        lbls[k].set_color(color)
    ax.set_ylim(-0.7, ymax - 0.3)
    ax.set_xlim(0, xmax * 1.15)
    ax.set_title(f"{prim}  (n={n_prot}, {round(100 * n_prot / N)}%)", color=color,
                 fontweight="bold", fontsize=10.5, loc="left")
    for s in ("top", "right", "left", "bottom"):   # no axes box; x-axis removed
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([])                              # bar-end labels carry the counts


def plot_tertiary_composition(path):
    prim_order = d[PRIM].value_counts().index.tolist()        # size desc
    grid = [prim_order[0:2], prim_order[2:4], prim_order[4:6]]  # 3 rows x 2 cols
    # each grid-row's panels share a y-scale = the taller panel's row count
    row_rows = [max(len(_panel_rows(p)[0]) for p in r) for r in grid]
    row_h = [max(h, 4) for h in row_rows]

    fig = plt.figure(figsize=(13.5, 0.34 * sum(row_h) + 1.6))
    gs = fig.add_gridspec(3, 2, height_ratios=row_h, hspace=0.22, wspace=0.55)
    for gi, gr in enumerate(grid):
        for gj, prim in enumerate(gr):
            _draw_panel(fig.add_subplot(gs[gi, gj]), prim, row_h[gi])
    fig.suptitle(f"syn3A proteome — tertiary function composition  (n = {N} proteins)",
                 fontsize=13, fontweight="bold")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


plot_tertiary_composition(OUT_FIG)


# ── Interactive HTML: click a tertiary function -> filter the protein table ───
def write_html(path):
    """Self-contained HTML (no external deps). Clickable composition bars on top;
    clicking a tertiary function filters the protein table below and scrolls to
    it. Table columns follow ORDER_COLS (incl. Protein Sequence). Headers sort on
    click; the current view downloads to CSV/Excel; the long sequence column is
    width-clamped (ellipsis, click a cell to expand) so it never blows out rows."""
    prim_order = d[PRIM].value_counts().index.tolist()
    cols = ORDER_COLS

    # records with NaN -> None so the JS table renders blanks
    records = json.loads(d[cols].astype(object).where(pd.notna(d[cols]), None)
                         .to_json(orient="records"))
    data_js = json.dumps(records).replace("</", "<\\/")
    cols_js = json.dumps(cols).replace("</", "<\\/")

    # build the clickable bar chart (same ordering as the PDF)
    panels = []
    for prim in prim_order:
        color = PRIM_COLORS.get(prim, "#777777")
        rows, n_prot = _panel_rows(prim)
        xmax = max([r[2] for r in rows if r[0] == "bar"], default=1)
        body, cur_sec = [], ""
        for r in rows:
            if r[0] == "sec":
                cur_sec = r[1]
                body.append(
                    f'<div class="sec-h clk" data-p="{esc(prim, quote=True)}" '
                    f'data-s="{esc(r[1], quote=True)}">{esc(r[1])} '
                    f'<span class="sn">({r[2]}, {round(100 * r[2] / N)}%)</span></div>')
            else:
                ter, c = r[1], r[2]
                w = 100 * c / xmax
                body.append(
                    f'<div class="bar" data-p="{esc(prim, quote=True)}" '
                    f'data-s="{esc(cur_sec, quote=True)}" data-t="{esc(ter, quote=True)}">'
                    f'<span class="lab">{esc(ter)}</span>'
                    f'<span class="track"><span class="fill" style="width:{w:.1f}%"></span></span>'
                    f'<span class="cnt">{c}</span></div>')
        panels.append(
            f'<section class="prim" style="--c:{color}">'
            f'<div class="prim-h clk" data-p="{esc(prim, quote=True)}">{esc(prim)} '
            f'<span class="pn">(n={n_prot}, {round(100 * n_prot / N)}%)</span></div>'
            f'{"".join(body)}</section>')
    chart_html = "".join(panels)

    html_doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>syn3A proteome — tertiary function composition</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
         margin: 0 auto; max-width: 1200px; padding: 20px; color:#222; }}
  h1 {{ font-size: 20px; }}
  .hint {{ color:#666; font-size:13px; margin-bottom:14px; }}
  .chart {{ display:grid; grid-template-columns: 1fr 1fr; gap: 10px 26px; }}
  @media (max-width: 800px) {{ .chart {{ grid-template-columns: 1fr; }} }}
  .prim-h {{ font-weight:700; color:var(--c); margin:10px 0 4px; font-size:15px; }}
  .pn {{ color:#888; font-weight:400; }}
  .sec-h {{ font-weight:700; color:var(--c); font-size:12px; margin:8px 0 2px 4px; opacity:.85; }}
  .sn {{ color:#999; font-weight:400; }}
  .clk {{ cursor:pointer; border-radius:4px; padding:1px 4px; }}
  .clk:hover {{ background:#f0f4fa; }}
  .bar {{ display:grid; grid-template-columns: 180px 1fr 30px; align-items:center;
          gap:6px; cursor:pointer; padding:1px 4px; border-radius:4px; }}
  .bar:hover {{ background:#f0f4fa; }}
  .bar.sel, .clk.sel {{ background:#fde9c8; }}
  .lab {{ font-size:11.5px; text-align:right; white-space:nowrap; overflow:hidden;
          text-overflow:ellipsis; }}
  .track {{ background:#eee; height:13px; border-radius:3px; }}
  .fill {{ display:block; height:100%; background:var(--c); border-radius:3px; }}
  .cnt {{ font-size:11px; color:#444; }}
  #panel {{ margin-top:26px; border-top:2px solid #ddd; padding-top:12px; }}
  #cap {{ font-weight:700; font-size:15px; }}
  #cap small {{ font-weight:400; color:#777; }}
  button {{ margin-left:10px; font-size:12px; padding:3px 8px; cursor:pointer; }}
  .tablebox {{ overflow-x:auto; margin-top:8px; max-height:75vh; }}
  table {{ border-collapse:collapse; font-size:12px; width:100%; }}
  th,td {{ border:1px solid #e2e2e2; padding:3px 7px; text-align:left;
           white-space:nowrap; vertical-align:top; }}
  th {{ position:sticky; top:0; background:#f5f5f5; cursor:pointer; user-select:none;
        white-space:nowrap; }}
  th.asc::after {{ content:" \\25B2"; font-size:9px; color:#888; }}
  th.desc::after {{ content:" \\25BC"; font-size:9px; color:#888; }}
  tr:nth-child(even) td {{ background:#fafafa; }}
  /* freeze the first three columns (Locus Tag, Gene Name, Gene Product) when
     the table is scrolled horizontally. Widths are fixed so the cumulative left
     offsets line up. Overflow truncates with an ellipsis; the title attribute
     on each cell reveals the full text on hover. */
  th:nth-child(1), td:nth-child(1) {{ position:sticky; left:0;
       min-width:140px; max-width:140px; overflow:hidden; text-overflow:ellipsis; }}
  th:nth-child(2), td:nth-child(2) {{ position:sticky; left:140px;
       min-width:90px;  max-width:90px;  overflow:hidden; text-overflow:ellipsis; }}
  th:nth-child(3), td:nth-child(3) {{ position:sticky; left:230px;
       min-width:320px; max-width:320px; overflow:hidden; text-overflow:ellipsis;
       border-right:2px solid #ccc; }}
  td:nth-child(1), td:nth-child(2), td:nth-child(3) {{ background:white; z-index:1; }}
  tr:nth-child(even) td:nth-child(1), tr:nth-child(even) td:nth-child(2),
  tr:nth-child(even) td:nth-child(3) {{ background:#fafafa; }}
  th:nth-child(1), th:nth-child(2), th:nth-child(3) {{ background:#f5f5f5; z-index:3; }}
  td.seq {{ font-family: ui-monospace, Menlo, Consolas, monospace; max-width:150px;
            overflow:hidden; text-overflow:ellipsis; white-space:nowrap; cursor:zoom-in; }}
  td.seq.exp {{ white-space:normal; word-break:break-all; max-width:340px; cursor:zoom-out; }}
</style></head>
<body>
<h1>syn3A proteome — tertiary function composition <small>(n = {N} proteins)</small></h1>
<div class="hint">Click any Primary / Secondary / Tertiary function to list its proteins below. Click a column header to sort; click a sequence cell to expand it.</div>
<div class="chart">{chart_html}</div>

<div id="panel">
  <span id="cap">Click a function to list its proteins</span>
  <button onclick="showAll()">Show all</button>
  <button onclick="downloadCSV()">Download CSV</button>
  <button onclick="downloadExcel()">Download Excel</button>
  <div class="tablebox"><table id="tbl"><thead><tr></tr></thead><tbody></tbody></table></div>
</div>

<script>
const DATA = {data_js};
const COLS = {cols_js};
const PRIM = "{esc(PRIM)}", SEC = "{esc(SEC)}", TER = "{esc(TER)}", SEQCOL = "{esc(SEQCOL)}";
const tb = document.querySelector("#tbl tbody");
const thr = document.querySelector("#tbl thead tr");
const cap = document.getElementById("cap");
let currentRows = DATA.slice(), curCaption = "All proteins";
let visibleCols = COLS.slice();
let sortCol = null, sortDir = 1;   // sortCol is a column NAME (robust to hiding)

function escHTML(v) {{ return (v===null||v===undefined) ? "" :
    String(v).replace(/&/g,"&amp;").replace(/</g,"&lt;"); }}
function cellHTML(col, v) {{
  const e = escHTML(v);
  const tt = e.replace(/"/g, '&quot;');         // hover tooltip = full value
  return col===SEQCOL ? '<td class="seq" title="'+tt+'">'+e+'</td>'
                      : '<td title="'+tt+'">'+e+'</td>';
}}
function draw() {{
  cap.innerHTML = curCaption + ' <small>(' + currentRows.length + ' proteins)</small>';
  thr.innerHTML = visibleCols.map(c => {{
    const cls = (c===sortCol) ? (sortDir>0 ? "asc" : "desc") : "";
    return '<th class="'+cls+'">'+escHTML(c)+'</th>';
  }}).join("");
  tb.innerHTML = currentRows.map(r =>
    "<tr>" + visibleCols.map(c => cellHTML(c, r[c])).join("") + "</tr>").join("");
}}
function applySort() {{
  if (sortCol !== null && visibleCols.indexOf(sortCol) !== -1) {{
    const c = sortCol;
    currentRows.sort((a, b) => {{
      let x = a[c], y = b[c];
      const xn = parseFloat(x), yn = parseFloat(y);
      const num = x!=null && y!=null && String(x).trim()!=="" && String(y).trim()!=="" &&
                  !isNaN(xn) && !isNaN(yn);
      const r = num ? (xn - yn)
                    : String(x==null?"":x).localeCompare(String(y==null?"":y));
      return r * sortDir;
    }});
  }}
  draw();
}}
// rows = filtered records; caption = HTML caption; hide = column names to omit
function setView(rows, caption, hide) {{
  currentRows = rows.slice(); curCaption = caption;
  visibleCols = COLS.filter(c => hide.indexOf(c) === -1);
  applySort();
}}
function clearSel() {{ document.querySelectorAll(".sel").forEach(e => e.classList.remove("sel")); }}
function jump() {{ document.getElementById("panel").scrollIntoView({{behavior:"smooth", block:"start"}}); }}
// reset to the default row order (the xlsx/genome order) by clearing any sort
function showAll() {{ clearSel(); sortCol = null; sortDir = 1; setView(DATA, "All proteins", []); }}

// click a Primary header -> filter to it, hide the Primary column
document.querySelectorAll(".prim-h.clk").forEach(el => el.addEventListener("click", () => {{
  clearSel(); el.classList.add("sel");
  const p = el.dataset.p;
  setView(DATA.filter(r => r[PRIM] === p), p, [PRIM]);
  jump();
}}));
// click a Secondary header -> filter, hide Primary + Secondary columns
document.querySelectorAll(".sec-h.clk").forEach(el => el.addEventListener("click", () => {{
  clearSel(); el.classList.add("sel");
  const p = el.dataset.p, s = el.dataset.s;
  setView(DATA.filter(r => r[PRIM] === p && r[SEC] === s), p + " &rsaquo; " + s, [PRIM, SEC]);
  jump();
}}));
// click a Tertiary bar -> filter, hide Primary + Secondary columns
document.querySelectorAll(".bar").forEach(bar => bar.addEventListener("click", () => {{
  clearSel(); bar.classList.add("sel");
  const p = bar.dataset.p, s = bar.dataset.s, t = bar.dataset.t;
  setView(DATA.filter(r => r[PRIM] === p && r[TER] === t && (!s || r[SEC] === s)),
          p + (s ? " &rsaquo; " + s : "") + " &rsaquo; " + t, [PRIM, SEC]);
  jump();
}}));
// sort on header click (column identified by name -> survives column hiding)
thr.addEventListener("click", e => {{
  const th = e.target.closest("th"); if (!th) return;
  const c = visibleCols[Array.prototype.indexOf.call(thr.children, th)];
  if (sortCol === c) {{ sortDir = -sortDir; }} else {{ sortCol = c; sortDir = 1; }}
  applySort();
}});
// expand a sequence cell on click
tb.addEventListener("click", e => {{
  if (e.target.classList.contains("seq")) e.target.classList.toggle("exp");
}});
// downloads of the current filtered+sorted rows (always the full column set)
function dl(blob, fname) {{
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = fname; a.click();
  URL.revokeObjectURL(a.href);
}}
function csvCell(v) {{ v = (v==null) ? "" : String(v);
  return /[",\\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v; }}
function downloadCSV() {{
  const lines = [COLS.map(csvCell).join(",")];
  currentRows.forEach(r => lines.push(COLS.map(c => csvCell(r[c])).join(",")));
  dl(new Blob(["\\ufeff" + lines.join("\\r\\n")], {{type:"text/csv;charset=utf-8"}}),
     "syn3A_proteins.csv");
}}
function downloadExcel() {{
  let h = "<table><tr>" + COLS.map(c => "<th>" + escHTML(c) + "</th>").join("") + "</tr>";
  currentRows.forEach(r => h += "<tr>" + COLS.map(c => "<td>" + escHTML(r[c]) + "</td>").join("") + "</tr>");
  h += "</table>";
  dl(new Blob(['\\ufeff<html><head><meta charset="utf-8"></head><body>' + h + "</body></html>"],
     {{type:"application/vnd.ms-excel"}}), "syn3A_proteins.xls");
}}
setView(DATA, "All proteins", []);
</script>
</body></html>"""
    with open(path, "w") as fh:
        fh.write(html_doc)


write_html(OUT_HTML)
say(f"\nSaved: {OUT_REPORT}")
# say(f"Saved: {OUT_TABLE}")
say(f"Saved: {OUT_XLSX}")
say(f"Saved: {OUT_FIG}")
say(f"Saved: {OUT_HTML}")
with open(OUT_REPORT, "w") as fh:
    fh.write("\n".join(_R) + "\n")
