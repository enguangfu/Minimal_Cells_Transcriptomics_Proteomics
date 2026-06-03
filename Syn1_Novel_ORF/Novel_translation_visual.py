import math
from typing import List, Tuple
import pandas as pd
import pandas as pd
import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


ISOFORM_TSV    = '../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv'
CONTEXT_PAD    = 2000   # flanking bp shown either side of ORF
MIN_ISO_READS  = 10      # drop isoforms below this read count
MAX_ISOFORMS   = 120    # cap per panel
ISO_LW_SCALE   = 0.5    # linewidth = scale * sqrt(n_reads)
ISO_GAP        = 1      # blank rows between TSS groups

DIGEST_TSV     = 'trypsin_digest/novel_peptide_digest_summary.tsv'   # regenerated on revised clusters
GFF_PATH       = '../Genomes_Input/syn1.genes.gff3'
MASS_SPEC_XLSX = 'Mass_Spec/Syn1.0_newSearch_20260402.xlsx'


def _load_genes_gff3(path):
    """Gene model with genomic coords + rna_type (replaces the now-dead
    syn1_genes_transcriptomics_proteomics.csv, which lacked genome coordinates)."""
    rows = []
    for line in open(path):
        if not line.strip() or line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) != 9 or p[2] != 'gene':
            continue
        attrs = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
        rows.append({'chrom': p[0], 'start0': int(p[3]) - 1, 'end0': int(p[4]), 'strand': p[6],
                     'locus_tag': attrs.get('locus_tag', ''),
                     'gene_name': (attrs.get('Name') or attrs.get('gene') or None),
                     'rna_type': attrs.get('rna_type', '')})
    return pd.DataFrame(rows)


digest = pd.read_csv(DIGEST_TSV, sep='\t')
genes  = _load_genes_gff3(GFF_PATH)
ms     = pd.read_excel(MASS_SPEC_XLSX, sheet_name='Proteins_all')
PG_COL = 'PG.NrOfStrippedSequencesIdentified (Experiment-wide)'

# ── Focal MS-confirmed novel ORFs ────────────────────────────────────────────
# IMPORTANT: the candidate set was REGENERATED on the revised (post Apr-22)
# isoform clusters, so the NOVEL_PEP_* IDs below are the NEW IDs. The Spectronaut
# MS re-search (Mass_Spec/Syn1.0_newSearch_20260402.xlsx) was run on the OLD
# clusters, so its protein accessions still carry the OLD IDs. Map NEW -> OLD:
#   NOVEL_PEP_002 (revised) == NOVEL_PEP_002 (old MS search)  : 118 aa intergenic ORF (0591/0592 gap)
#   NOVEL_PEP_043 (revised) == NOVEL_PEP_030 (old MS search)  : 225 aa N-term extension of mmyCIVR/0768
# i.e. what is NOVEL_PEP_043 here USED TO BE NOVEL_PEP_030 in the proteome
# new-search Excel (old clusters). Old NOVEL_PEP_032 was a near-duplicate of 002
# and is dropped.
FOCUS_MS_ID = {'NOVEL_PEP_002': 'NOVEL_PEP_002',
               'NOVEL_PEP_043': 'NOVEL_PEP_030'}
FOCUS = list(FOCUS_MS_ID)
foc   = digest[digest['novel_peptide_id'].isin(FOCUS)].set_index('novel_peptide_id')


isoforms_all = pd.read_csv(ISOFORM_TSV, sep='\t')

# ── coordinate transform (mirrors OperonCoord.tx_of_genome_pos0) ─────────────
# Maps a genomic pos0 to a tx-space x where 5' end = small x, 3' end = large x.
# + strand: x = pos0 - win_s   (low coord is 5')
# - strand: x = win_e - pos0   (high coord is 5')
def tx(pos0: int, win_s: int, win_e: int, strand: str) -> int:
    if strand == '+':
        return pos0 - win_s
    else:
        return win_e - pos0

def tx_to_genome(x: int, win_s: int, win_e: int, strand: str) -> int:
    if strand == '+':
        return win_s + x
    else:
        return win_e - x

# ── interval packing (identical to Operon_Visualization) ─────────────────────
def _pack_rows(intervals: List[Tuple[int, int]]) -> List[int]:
    """Greedy non-overlapping row assignment, sorted by left endpoint."""
    order = sorted(range(len(intervals)), key=lambda i: (intervals[i][0], intervals[i][1]))
    row_ends: List[int] = []
    rows_out = [0] * len(intervals)
    for i in order:
        a, b = intervals[i]
        placed = False
        for r, end in enumerate(row_ends):
            if end <= a:
                row_ends[r] = b
                rows_out[i] = r
                placed = True
                break
        if not placed:
            row_ends.append(b)
            rows_out[i] = len(row_ends) - 1
    return rows_out


def layout_isoforms(iso_df: pd.DataFrame,
                    win_s: int, win_e: int,
                    orf_s: int, orf_e: int,
                    strand: str) -> pd.DataFrame:
    """
    Adapted directly from layout_isoform_tracks in Operon_Visualization.

    Key conventions (matching that notebook):
      start_pos0 = pos5p0  (5' end; high coord on minus strand)
      end_pos0   = pos3p0  (3' end; low coord on minus strand)

    Filter: keep isoforms that fully cover [orf_s, orf_e]
      i.e. start0 <= orf_s  AND  end0 >= orf_e  (genomic low/high coords)

    Coordinate transform: tx() maps genomic pos -> tx-space where
      x=0 is the left edge of the window and x increases left->right,
      with 5' end always mapping to a SMALLER x than the 3' end.
      This is identical to oc.tx_of_genome_pos0 in Operon_Visualization.

    Grouping: by start_pos0 (= pos5p0), sort=False (preserves pre-sort order).
    """
    iso = iso_df.copy()
    
    # Use the same column names as Operon_Visualization for clarity
    iso['start_pos0'] = iso['pos5p0']   # 5' end
    iso['end_pos0']   = iso['pos3p0']   # 3' end

    # Sort exactly as Operon_Visualization does:
    #   + strand: start_pos0 asc, end_pos0 asc, n_reads desc
    #   - strand: start_pos0 desc, end_pos0 desc, n_reads desc
    if strand == '+':
        iso = iso.sort_values(['start_pos0', 'end_pos0', 'n_reads'],
                              ascending=[True, True, False])
    else:
        iso = iso.sort_values(['start_pos0', 'end_pos0', 'n_reads'],
                              ascending=[False, False, False])
    iso = iso.head(MAX_ISOFORMS).copy()

    # Convert to tx-space, clip to window, then take min/max for left/right
    tx_lefts, tx_rights = [], []
    for _, r in iso.iterrows():
        s = int(r['start_pos0'])
        e = int(r['end_pos0'])
        s_clip = min(max(s, win_s), win_e)
        e_clip = min(max(e, win_s), win_e)
        x0 = tx(s_clip, win_s, win_e, strand)
        x1 = tx(e_clip, win_s, win_e, strand)
        tx_lefts.append(min(x0, x1))
        tx_rights.append(max(x0, x1))
    iso['tx_left']  = tx_lefts
    iso['tx_right'] = tx_rights
    iso = iso[iso['tx_right'] > iso['tx_left']].copy()
    if iso.empty:
        return iso

    # Thickness + alpha (same as Operon_Visualization)
    def thickness(n):
        return ISO_LW_SCALE * math.sqrt(max(0.0, float(n)))

    def _scale_alpha(n, nmin, nmax):
        if nmax <= nmin:
            return 0.9
        return float(0.5 + 0.45 * (n - nmin) / (nmax - nmin))

    iso['y']     = float('nan')
    iso['lw']    = float('nan')
    iso['alpha'] = float('nan')

    base_y = 0
    # sort=False preserves the pre-sort order within each group
    for g_start, g in iso.groupby('start_pos0', sort=False):
        g = g.sort_values(['n_reads', 'tx_right'], ascending=[False, True]).copy()
        intervals = list(zip(g['tx_left'].astype(int), g['tx_right'].astype(int)))
        rows  = _pack_rows(intervals)
        nvals = g['n_reads'].astype(float).to_numpy()
        nmin, nmax = float(nvals.min()), float(nvals.max())
        iso.loc[g.index, 'y']     = [base_y + r for r in rows]
        iso.loc[g.index, 'lw']    = [max(0.3, thickness(float(n))) for n in nvals]
        iso.loc[g.index, 'alpha'] = [_scale_alpha(float(n), nmin, nmax) for n in nvals]
        base_y += (max(rows) + 1) + ISO_GAP

    iso['group_id'] = iso['start_pos0'].astype(int)
    return iso


def draw_isoforms(ax, iso_df: pd.DataFrame,
                  win_s: int, win_e: int,
                  strand: str, orf_s: int, orf_e: int):
    """
    Draw isoforms as FancyArrowPatch, pointing in transcript direction.
    x-axis is tx-space (0 = left edge of window, increases right).
    Arrow always points right (5'->3' in tx-space), matching Operon_Visualization.
    """
    iso = layout_isoforms(iso_df, win_s, win_e, orf_s, orf_e, strand)
    # print("After layout_isoforms, the iso df is")
    # print(iso[["isoform_id", "start_pos0", "end_pos0", "tx_left", "tx_right", "y", "lw", "alpha"]].head(10)) 
    if iso.empty:
        ax.text(0.5, 0.5, 'No isoforms', transform=ax.transAxes,
                ha='center', va='center', fontsize=9, color='gray')
        ax.set_yticks([])
        return

    starts = sorted(iso['group_id'].unique())
    cmap   = plt.get_cmap('tab10')
    cdict  = {s: cmap(i % cmap.N) for i, s in enumerate(starts)}
    win_len = win_e - win_s

    for _, r in iso.iterrows():
        xl    = float(r['tx_left'])
        xr    = float(r['tx_right'])
        y     = float(r['y'])
        col   = cdict[int(r['group_id'])]
        lw    = float(r['lw'])
        alpha = float(r['alpha'])
        width = max(1.0, xr - xl)

        # Arrow always points right in tx-space (5'->3'), same as Operon_Visualization
        ax.add_patch(mpatches.FancyArrowPatch(
            (xl, y), (xr, y),
            arrowstyle='-|>',
            linewidth=lw, color=col, alpha=alpha,
            shrinkA=0, shrinkB=0,
        ))

    y_max = float(iso['y'].max())
    ax.set_xlim(0, win_len)
    ax.set_ylim(-1, y_max + 2)
    ax.set_yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('RNA isoforms\n(grouped by TSS)', fontsize=8)

    # x-axis labels as genomic coordinates
    _ws2, _we2, _st2 = win_s, win_e, strand
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _, ws=_ws2, we=_we2, st=_st2: f'{tx_to_genome(int(x), ws, we, st)/1000:.2f} kb'))


# ── gene-overlap text report ─────────────────────────────────────────────────
def gene_overlaps(genes_df, chrom, start0, end0, orf_strand):
    mask = ((genes_df['chrom'] == chrom) & (genes_df['start0'] < end0) &
            (genes_df['end0'] > start0))
    hits = genes_df[mask].copy()
    hits['relationship'] = hits['strand'].apply(
        lambda s: 'sense' if s == orf_strand else 'antisense')
    return hits

print('=' * 72)
print('2. Translated region projected onto genome — gene overlaps')
print('=' * 72)
for pid in FOCUS:
    row    = foc.loc[pid]
    gs, ge = int(row['orf_genomic_start0']), int(row['orf_genomic_end0'])
    strand = row['strand']
    ms_id  = FOCUS_MS_ID[pid]   # MS Excel uses the OLD-cluster ID
    ms_row = ms[ms['PG.ProteinAccessions'].str.contains(ms_id, na=False)]
    pg_n   = ms_row[PG_COL].values[0] if len(ms_row) else 'n/d'
    id_note = '' if ms_id == pid else f'   [== {ms_id} in the old-cluster MS search Excel]'
    print(f'\n  {pid}  ({row["orf_aa_len"]} aa, {strand} strand){id_note}')
    print(f'    ORF coords: {gs:,}-{ge:,}  |  MS peptides (exp-wide, searched as {ms_id}): {pg_n}')
    overlap = gene_overlaps(genes, 'CP002027.1', gs, ge, strand)
    if overlap.empty:
        print('    Gene overlap: none (fully intergenic)')
    else:
        for _, gr in overlap.iterrows():
            ov_bp = min(ge, int(gr['end0'])) - max(gs, int(gr['start0']))
            name  = gr['gene_name'] if pd.notna(gr['gene_name']) else gr['locus_tag']
            print(f'    [{gr["relationship"].upper()}] {name} ({gr["locus_tag"]}, '
                  f'{gr["rna_type"]})  {int(gr["start0"]):,}-{int(gr["end0"]):,}  '
                  f'{ov_bp} bp  ({ov_bp/(ge-gs)*100:.0f}% of ORF / '
                  f'{ov_bp/(gr["end0"]-gr["start0"])*100:.0f}% of gene)')

# ── figure ───────────────────────────────────────────────────────────────────
GENE_COLORS = {'mRNA': '#4393C3', 'pseudo': '#B2ABD2', 'rRNA': '#E08214',
               'tRNA': '#7FBC41', 'tmRNA': '#E08214', 'ncRNA': '#DFC27D'}
GENE_COLORS_STRAND  = {'+': '#1f77b4', '-': '#ff7f0e'}  # blue/orange by strand
ORF_COLORS  = {'NOVEL_PEP_002': '#D6604D', 'NOVEL_PEP_043': '#A50026'}
GENE_LABEL_FONTSIZE = 6.5


def draw_gene_arrows(ax, genes_df: pd.DataFrame,
                     win_s: int, win_e: int, strand: str):
    """
    Draw genes as polygon arrows (rectangle + triangular head), mirroring
    draw_gene_arrows in Operon_Visualization.py.

    X-axis is tx-space (0 = window left, increases rightward, 5'→3').
    Y_BASE = 1.0.  Ylim, yticks, and spines are left to the caller so the
    ORF row drawn at y=2 on the same axis remains visible.
    """
    Y_BASE    = 1.0
    RECT_H    = 0.28
    TRI_H     = 0.28
    HEAD_FRAC = 0.18
    HEAD_MIN  = 20

    xlim_lo, xlim_hi = ax.get_xlim()
    ax.hlines(Y_BASE, xlim_lo, xlim_hi, color='black', lw=2, zorder=1)

    if genes_df.empty:
        ax.text(0.01, Y_BASE, 'No genes in interval', va='center', fontsize=7)
        return

    for _, r in genes_df.sort_values('start0').reset_index(drop=True).iterrows():
        g0      = int(r['start0'])
        g1      = int(r['end0'])
        gstrand = str(r['strand'])

        x0     = tx(g0, win_s, win_e, strand)
        x1     = tx(g1, win_s, win_e, strand)
        xleft  = min(x0, x1)
        xright = max(x0, x1)
        xcen   = (xleft + xright) / 2
        width  = xright - xleft

        head_len = min(max(HEAD_MIN, width * HEAD_FRAC), width)

        color        = GENE_COLORS_STRAND.get(gstrand, '#AAAAAA')
        points_right = (gstrand == strand)

        if points_right:
            tip_x  = xright
            base_x = xright - head_len
            verts = [
                (xleft,  Y_BASE - RECT_H / 2),
                (base_x, Y_BASE - RECT_H / 2),
                (base_x, Y_BASE - TRI_H  / 2),
                (tip_x,  Y_BASE),
                (base_x, Y_BASE + TRI_H  / 2),
                (base_x, Y_BASE + RECT_H / 2),
                (xleft,  Y_BASE + RECT_H / 2),
            ]
        else:
            tip_x  = xleft
            base_x = xleft + head_len
            verts = [
                (xright, Y_BASE - RECT_H / 2),
                (base_x, Y_BASE - RECT_H / 2),
                (base_x, Y_BASE - TRI_H  / 2),
                (tip_x,  Y_BASE),
                (base_x, Y_BASE + TRI_H  / 2),
                (base_x, Y_BASE + RECT_H / 2),
                (xright, Y_BASE + RECT_H / 2),
            ]

        ax.add_patch(mpatches.Polygon(
            verts, closed=True,
            facecolor=color, edgecolor='black',
            lw=0, alpha=0.85, zorder=2,
        ))

        gene_name = r.get('gene_name', '')
        locus_tag = str(r.get('locus_tag', ''))
        label = str(gene_name) if pd.notna(gene_name) and str(gene_name).strip() else locus_tag
        if label:
            ax.text(xcen, Y_BASE + TRI_H / 2 + 0.06, label,
                    ha='center', va='bottom',
                    fontsize=GENE_LABEL_FONTSIZE, color=color, clip_on=True)

for i, pid in enumerate(FOCUS):
    fig, axes = plt.subplots(2, 1, figsize=(15, (1 + 1) * 2.5),
                             gridspec_kw={'height_ratios': [1, 1]},
                             sharex=False)
    ax_gene = axes[0]
    ax_iso  = axes[1]

    row    = foc.loc[pid]
    gs, ge = int(row['orf_genomic_start0']), int(row['orf_genomic_end0'])
    strand = row['strand']
    chrom  = 'CP002027.1'
    win_s  = gs - CONTEXT_PAD
    win_e  = ge + CONTEXT_PAD
    win_len = win_e - win_s
    print("gs, ge:", gs, ge)
    print("win_s, win_e:", win_s, win_e)
    # ── gene track (tx-space x-axis) ─────────────────────────────────────
    ctx = genes[
        (genes['chrom']  == chrom) &
        (genes['strand'] == strand) &
        (genes['start0']  > win_s) &
        (genes['end0']   < win_e)
    ]
    ax_gene.set_xlim(0, win_len)
    ax_gene.set_ylim(0, 2.8)

    GENE_Y = 1.0
    draw_gene_arrows(ax_gene, ctx, win_s, win_e, strand)

    ORF_Y = 2;  ORF_H = 0.4
    orf_color = ORF_COLORS.get(pid, '#E04444')
    # ORF coords in tx-space
    ox0 = tx(gs, win_s, win_e, strand)
    ox1 = tx(ge, win_s, win_e, strand)
    oxl, oxr = min(ox0, ox1), max(ox0, ox1)
    ax_gene.add_patch(mpatches.FancyBboxPatch(
        (oxl, ORF_Y - ORF_H/2), oxr - oxl, ORF_H,
        boxstyle='round,pad=0', lw=1.2,
        edgecolor='#800000', facecolor=orf_color, alpha=0.9, zorder=3))
    ax_gene.text((oxl + oxr) / 2, ORF_Y + ORF_H/2 + 0.06,
                 f'{pid}  ({row["orf_aa_len"]} aa)',
                 ha='center', va='bottom', fontsize=8, fontweight='bold',
                 color='#800000', clip_on=True)
    for xv in [oxl, oxr]:
        ax_gene.axvline(xv, color='#800000', lw=0.9, ls='--', alpha=0.5, zorder=4)

    ax_gene.set_yticks([GENE_Y, ORF_Y])
    ax_gene.set_yticklabels([f'genes ({strand})', 'Novel ORF'], fontsize=7)
    _ws, _we, _st = win_s, win_e, strand
    ax_gene.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _, ws=_ws, we=_we, st=_st: f'{tx_to_genome(int(x), ws, we, st)/1000:.2f} kb'))
    ax_gene.tick_params(axis='x', labelbottom=False)
    ms_id   = FOCUS_MS_ID[pid]
    id_note = '' if ms_id == pid else f' (MS search: {ms_id}, old clusters)'
    ax_gene.set_title(
        f'{pid}{id_note}  [{strand} strand, ORF {gs:,}–{ge:,}, {row["orf_aa_len"]} aa]',
        fontsize=9, fontweight='bold')
    gene_color = GENE_COLORS_STRAND.get(strand, '#AAAAAA')
    ax_gene.legend(handles=[
        mpatches.Patch(color=gene_color, label=f'gene ({strand} strand)'),
        mpatches.Patch(color=orf_color,  label='Novel ORF'),
    ], fontsize=7, loc='upper right', framealpha=0.85)

    # ── isoform track ────────────────────────────────────────────────────
    iso_win = isoforms_all[
        (isoforms_all['chrom']  == chrom) &
        (isoforms_all['strand'] == strand) &
        (isoforms_all['start0'] <= gs) &
        (isoforms_all['n_reads'] >= MIN_ISO_READS) &
        (isoforms_all['end0']   >= ge)
    ].copy()
    
    print(f'Isoforms fully covering ORF with ≥{MIN_ISO_READS} reads: {len(iso_win)}')
    
    draw_isoforms(ax_iso, iso_win, win_s, win_e, strand, gs, ge)

    for xv in [oxl, oxr]:
        ax_iso.axvline(xv, color='#800000', lw=0.9, ls='--', alpha=0.5, zorder=5)

    ax_iso.set_title(
        f'{iso_win.shape[0]} isoforms (≥{MIN_ISO_READS} reads) fully covering ORF  '
        f'[coloured by TSS / pos5p0]',
        fontsize=8)

    plt.tight_layout()
    out_pdf = f'novel_pep_{pid}_genomic_context.pdf'
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    plt.show()
    print(f'Figure saved: {out_pdf}')
