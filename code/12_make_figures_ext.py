#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 图 5-7（扩展图组；图内全英文——规则 152）。

Fig5 LOQO stability: Tier A probe survival rate per dropped queue (from loqo_stability.json)
Fig6 Contrast layer: competitor gene-based hits vs our high-confidence layer (intersection 1)
Fig7 High-confidence lead SNP chromosome distribution (BTA 1-29)
高置信层规模由 integrated 表实时统计（strict 口径 = PIP>=0.5 & p<0.05），不硬编码。
"""
import json, os, sys, csv
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, 'output')
FIG = os.path.join(HERE, 'figures')
plt.rcParams.update({'font.size': 10, 'axes.linewidth': 0.8})


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, f'{name}.{ext}'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('saved', name)


def fig5_loqo():
    loqo = json.load(open(os.path.join(OUT, 'loqo_stability.json')))
    qs = list(loqo)
    rates = [loqo[q]['rate'] or 0 for q in qs]
    labels = [q.replace('GSE', '') for q in qs]
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    bars = ax.bar(labels, rates, color=['#2ca02c' if r >= 0.9 else '#f0ad4e' if r >= 0.75 else '#d9534f'
                                        for r in rates], edgecolor='black', lw=0.5)
    for b, q in zip(bars, qs):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.01,
                f"{loqo[q]['survived']}/{loqo[q]['tested']}", ha='center', fontsize=7.5)
    ax.axhline(1.0, color='gray', ls=':', lw=0.8)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel('Tier A survival rate after dropping one queue')
    ax.set_xlabel('Left-out queue')
    ax.set_title('Leave-one-queue-out stability of fully direction-consistent signals',
                 fontsize=11, fontweight='bold')
    save(fig, 'Fig3_loqo_stability')


def load_hi():
    """高置信层（PIP>=0.5 & p<0.05）实时统计，避免图内硬编码数字漂移。"""
    rows = [r for r in csv.DictReader(open(os.path.join(OUT, 'de_eqtl_pip_integrated.tsv'),
                                           encoding='utf-8'), delimiter='\t')
            if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05]
    return rows


def fig6_contrast():
    bl = json.load(open(os.path.join(OUT, 'mini_qtl_gene_baseline.json')))
    rec_path = os.path.join(OUT, 'competitor_reconciled.json')
    if os.path.exists(rec_path):
        rec = json.load(open(rec_path))
        n_comp = rec['unique_reconciled_ensembl']
        it = set(json.load(open(os.path.join(OUT, 'mini_qtl_intersection_reconciled.json'))))
    else:
        n_comp = len(bl)
        it = set(json.load(open(os.path.join(OUT, 'mini_qtl_intersection.json'))))
    n_hi = len(load_hi())
    fig, ax = plt.subplots(figsize=(7.0, 4.9))
    bars = ax.bar([0, 1, 2], [n_comp, n_hi, len(it)],
                  color=['#4a90d9', '#2ca02c', '#d9534f'], width=0.58,
                  edgecolor='black', linewidth=0.5)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([
        f'Competitor gene-based\n({len(bl)} identifiers ->\n{n_comp} Ensembl genes)',
        f'This study\n(high-confidence tier,\n{n_hi} genes)',
        f'Shared\n({len(it)} genes)'], fontsize=9)
    for x, v in zip((0, 1, 2), (n_comp, n_hi, len(it))):
        ax.text(x, v + n_hi * 0.02, str(v), ha='center', va='bottom',
                fontweight='bold', fontsize=11)
    ax.set_ylim(0, n_hi * 1.18)
    ax.set_ylabel('Number of genes', fontsize=10)
    ax.tick_params(axis='y', labelsize=9)
    ax.grid(axis='y', color='#dddddd', lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title('Two evidence streams are near-disjoint', fontsize=11, fontweight='bold', pad=10)
    save(fig, 'Fig7_contrast_layers')


def fig7_chrom():
    hi = load_hi()
    snps = [r['lead_SNP'] for r in hi]
    chrn = Counter(int(s.split('_')[0]) for s in snps)
    xs = list(range(1, 30))
    ys = [chrn.get(x, 0) for x in xs]
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    colors = ['#4a90d9' if i % 2 == 0 else '#1f5c99' for i in range(29)]
    ax.bar(xs, ys, color=colors, edgecolor='black', lw=0.4)
    for x, y in zip(xs, ys):
        if y:
            ax.text(x, y + 0.3, str(y), ha='center', fontsize=6.5)
    ax.set_xticks(xs)
    ax.set_xlabel('Bos taurus autosome (BTA)')
    ax.set_ylabel('High-confidence lead SNPs')
    ax.set_title(f'Genomic distribution of {len(hi)} high-confidence lead SNPs (max PIP >= 0.5, p < 0.05)',
                 fontsize=11, fontweight='bold')
    save(fig, 'Fig5_lead_snp_chromosomes')


if __name__ == '__main__':
    fig5_loqo()
    fig6_contrast()
    fig7_chrom()
