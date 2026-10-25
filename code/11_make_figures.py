#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 图 1-4 生成（matplotlib 统一管线 300dpi PNG+PDF；图内文本全英文——规则 152）。

数字取 2026-10-07 定稿口径（strict：整合层 = DE∩eQTL显著∩finemapping）：
  60 GSE -> 9 队列（339 samples）-> 7 芯片+2 RNA-seq -> 24,128 探针 -> 1,886/24,597 基因
  -> 3,395 整合基因 / 286 高置信层
Fig2  合并基因级 meta 火山图（beta vs p_DL，标注三源基因与 BH 显著；Bonferroni 参考线）
Fig3  DE-eQTL 整合散点（max_PIP vs -log10 p_DL；286 高置信层高亮）
Fig4  Top 候选四源打分条形图（Top 20 堆叠 S1/S3/S4 贡献）
"""
import csv, json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, 'output')
FIG = os.path.join(HERE, 'figures')
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({'font.size': 10, 'axes.linewidth': 0.8})


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, f'{name}.{ext}'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('saved', name)


def fig1_data_flow():
    stages = [
        ('GEO mastitis pool screened\n(60 expression series)', 60),
        ('mRNA DE cohorts assembled\n(9 series, 339 samples)', 9),
        ('Quantitative synthesis\n(7 cohorts, 257 samples: 5 chip + 2 RNA-seq)', 7),
        ('Cross-queue probe meta\n(24,128 probes; five GPL2112 cohorts)', 24128),
        ('Gene-level combined meta\n(24,597 Ensembl genes)', 24597),
        ('Integrated with mammary eQTL + PIP\n(3,395 genes; 286 high-confidence)', 3395),
    ]
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    n = len(stages)
    for i, (label, v) in enumerate(stages):
        y = n - 1 - i
        frac = 1.0 if i < 3 else (v / 24597) ** 0.3
        width = 0.78 * frac
        left = 0.5 - width / 2
        ax.barh(y, width, left=left, height=0.6,
                color=plt.cm.Blues(0.25 + 0.62 * i / (n - 1)), edgecolor='black', lw=0.6)
        ax.text(0.5, y, f'{v:,}', ha='center', va='center', fontsize=11,
                fontweight='bold', color='white' if i >= 3 else 'black')
        ax.text(-0.02, y, label, ha='right', va='center', fontsize=9)
        if i:
            ax.annotate('', xy=(0.5, y + 0.42), xytext=(0.5, y + 0.78),
                        arrowprops=dict(arrowstyle='-|>', lw=1.1, color='gray'))
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.6, n - 0.3)
    ax.axis('off')
    ax.set_title('Analysis data flow: GEO queues to four-source prioritization',
                 fontsize=12, fontweight='bold', pad=12)
    save(fig, 'Fig1_data_flow')


def fig2_volcano():
    rows = list(csv.DictReader(open(os.path.join(OUT, 'de_meta_gene_combined.tsv'), encoding='utf-8'), delimiter='\t'))
    beta = np.array([float(r['beta_meta']) for r in rows])
    p = np.array([float(r['p_DL']) for r in rows])
    src = [r['sources'] for r in rows]
    # BH over genes with a defined combined p only (21,442), matching the manuscript
    ok = np.isfinite(p)
    q = p[ok]
    order = np.argsort(q)
    bh_ok = np.empty_like(q)
    run = 1.0
    for rank, i in enumerate(order[::-1], 1):
        run = min(run, q[i] * len(q) / (len(q) - rank + 1))
        bh_ok[i] = run
    bh = np.full_like(p, np.nan)
    bh[ok] = bh_ok
    x = beta
    y = -np.log10(np.clip(p, 1e-300, 1))
    multi = np.array(['chip' in s and 'GSE13' in s and 'GSE75' in s for s in src])
    hi = np.isfinite(bh) & (bh < 0.05)
    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    ax.scatter(x[~hi], y[~hi], s=4, c='#c9c9c9', alpha=0.5, lw=0)
    ax.scatter(x[hi], y[hi], s=6, c='#4a90d9', alpha=0.55, lw=0,
               label=f'BH q<0.05 ({int(hi.sum())})')
    ax.scatter(x[multi], y[multi], s=10, c='#d9544a', alpha=0.8, lw=0,
               label=f'three-source genes ({int(multi.sum())})')
    named = {'ENSBTAG00000011586': 'RNASE1', 'ENSBTAG00000020313': 'PRL',
             'ENSBTAG00000005660': 'FOS', 'ENSBTAG00000022183': 'NFKBIA',
             'ENSBTAG00000008102': 'MET', 'ENSBTAG00000008463': 'CD3G'}
    for r, xi, yi in zip(rows, x, y):
        if r['gene'] in named and r['gene'] in set(rows[i]['gene'] for i in np.where(multi)[0]):
            ax.annotate(named[r['gene']], (xi, yi), fontsize=7.5, ha='center', va='bottom')
    ax.set_xlabel('Combined meta beta (infection effect)')
    ax.set_ylabel('-log10(p) random-effects meta')
    ax.axhline(-np.log10(0.05 / int(ok.sum())), color='gray', ls=':', lw=0.8)
    ax.text(ax.get_xlim()[0] + 0.1, -np.log10(0.05 / int(ok.sum())) + 0.15, 'Bonferroni', fontsize=7, color='gray')
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    ax.set_title('Gene-level combined meta-analysis (seven synthesized cohorts)', fontsize=11, fontweight='bold')
    save(fig, 'Fig2_volcano_combined_meta')


def fig3_pip():
    rows = list(csv.DictReader(open(os.path.join(OUT, 'de_eqtl_pip_integrated.tsv'), encoding='utf-8'), delimiter='\t'))
    pip = np.array([float(r['max_PIP']) for r in rows])
    p = np.array([float(r['p_DL']) for r in rows])
    hi = (pip >= 0.5) & (p < 0.05)
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    ax.scatter(pip[~hi], -np.log10(np.clip(p[~hi], 1e-300, 1)), s=5, c='#c9c9c9', alpha=0.5, lw=0)
    ax.scatter(pip[hi], -np.log10(np.clip(p[hi], 1e-300, 1)), s=9, c='#2ca02c', alpha=0.7, lw=0,
               label=f'high-confidence (PIP>=0.5 & p<0.05): {int(hi.sum())}')
    ax.axvline(0.5, color='gray', ls='--', lw=0.8)
    ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
    ax.set_xlabel('Maximum finemapping PIP (mammary cis-eQTL)')
    ax.set_ylabel('-log10(p) combined DE meta')
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    ax.set_title('DE evidence x mammary eQTL finemapping confidence', fontsize=11, fontweight='bold')
    save(fig, 'Fig4_DE_eQTL_PIP')


def fig4_scoring():
    rows = list(csv.DictReader(open(os.path.join(OUT, 'top_candidates.tsv'), encoding='utf-8'), delimiter='\t'))[:20]
    respath = os.path.join(OUT, 'ensembl_symbol_resolution.json')
    symres = json.load(open(respath, encoding='utf-8'))['mapping'] if os.path.exists(respath) else {}
    labels = [(r['symbol'] or symres.get(r['gene'], {}).get('symbol')
               or r['gene'].replace('ENSBTAG000000', 'ENS')) for r in rows]
    s1 = np.array([float(r['S1_DE']) for r in rows])
    s3 = np.array([float(r['S3_PIP']) for r in rows])
    s4 = np.array([float(r['S4_network']) for r in rows])
    y = np.arange(len(rows))[::-1]
    fig, ax = plt.subplots(figsize=(7.8, 6.8))
    h1 = ax.barh(y, s1 * 0.35, height=0.72, color='#4a90d9', label='DE meta (w=0.35)')
    h3 = ax.barh(y, s3 * 0.30, height=0.72, left=s1 * 0.35, color='#2ca02c',
                 label='fine-mapping PIP (w=0.30)')
    from matplotlib.patches import Patch
    handles = [Patch(facecolor='#4a90d9', label='DE meta (w=0.35)'),
               Patch(facecolor='#2ca02c', label='fine-mapping PIP (w=0.30)')]
    if (s4 > 0).any():
        ax.barh(y, s4 * 0.10, height=0.72, left=s1 * 0.35 + s3 * 0.30, color='#9467bd')
        handles.append(Patch(facecolor='#9467bd', label='STRING network (w=0.10)'))
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlim(0, max(0.6, (s1 * 0.35 + s3 * 0.30 + s4 * 0.10).max() * 1.08))
    ax.set_xlabel('Priority score (0-1; the eQTL-presence term of 0.25 is constant and not shown)',
                  fontsize=9.5)
    ax.tick_params(axis='x', labelsize=9)
    ax.grid(axis='x', color='#dddddd', lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(handles=handles, frameon=False, fontsize=8.5, ncol=len(handles),
              loc='upper center', bbox_to_anchor=(0.5, -0.11), handlelength=1.6,
              columnspacing=1.6)
    ax.set_title('Top 20 candidate genes by four-source priority scoring', fontsize=11,
                 fontweight='bold', pad=10)
    save(fig, 'Fig6_priority_scoring')


if __name__ == '__main__':
    fig1_data_flow()
    fig2_volcano()
    fig3_pip()
    fig4_scoring()
