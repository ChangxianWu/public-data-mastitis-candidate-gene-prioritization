#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 4b/4c: 四源加权打分表 + Top-N 候选基因决策表。

源与权重（0-1 归一后加权；QTL 层缺位以占位列保留，恢复后并入）：
  S1 DE 证据: -log10(p_DL_combined) 归一（合并 meta）
  S2 eQTL 存在: ∈ mammary significant (7,206) = 0.5 底分（整合层纳入判据，层内为常数）
  S3 精细定位: max_PIP
  S4 网络中心性: STRING 子网 weighted degree 归一
打分 = 0.35*S1 + 0.25*S2 + 0.30*S3 + 0.10*S4
高置信候选域：合并 meta p_DL<0.05 ∧ PIP>=0.5
输出 work/output/priority_scoring.tsv + top_candidates.tsv（Top 30 决策表）

2026-10-25（评审 M2/M6）：
  - 输出同时保留未加权度 STRING_degree 与加权度 STRING_weighted_degree（旧版把截断加权度
    标为 STRING_degree，与正文引用的未加权度不一致）
  - 验证建议按可执行性分层：无官方符号者不给 qPCR 引物建议
"""
import csv, io, json, math, os, sys
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.normpath(os.path.join(HERE, 'output'))


def outpath(name):
    base = os.path.normpath(OUT)
    p = os.path.normpath(base + os.sep + name)
    if not p.startswith(base + os.sep):
        raise SystemExit('output path escaped OUT dir: ' + name)
    return p


def write_tsv(name, fieldnames, rows):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
    w.writeheader()
    for r in rows:
        w.writerow(r)
    Path(outpath(name)).write_text(buf.getvalue(), encoding='utf-8')


def write_json(name, obj):
    Path(outpath(name)).write_text(json.dumps(obj, indent=1), encoding='utf-8')


def main():
    base = {}
    n_nan = 0
    for r in csv.DictReader(open(outpath('de_eqtl_pip_integrated.tsv'), encoding='utf-8'), delimiter='\t'):
        g = r['gene']
        p = float(r['p_DL'])
        if not math.isfinite(p):
            n_nan += 1
            continue
        base[g] = {'symbol': r['symbol'], 'p_DL': p, 'pip': float(r['max_PIP']),
                   'n_queues': int(r['n_queues']), 'beta': float(r['beta_meta']), 'snp': r['lead_SNP']}
    print(f'integrated genes: {len(base) + n_nan}  ranked (finite meta-p): {len(base)}  excluded (undefined p): {n_nan}')
    wd, deg = {}, {}
    for r in csv.DictReader(open(outpath('string_degree.tsv'), encoding='utf-8'), delimiter='\t'):
        wd[r['gene']] = float(r['weighted_degree'])
        deg[r['gene']] = int(r['degree'])
    ps = [v['p_DL'] for v in base.values()]
    wds = list(wd.values())

    def norm(x, lo, hi):
        return (x - lo) / (hi - lo) if hi > lo else 0.0

    pmin, pmax = min(ps), max(ps)
    wmin, wmax = min(wds), max(wds)
    score_rows = []
    for g, v in base.items():
        s1 = norm(-math.log10(max(v['p_DL'], 1e-300)), -math.log10(pmax), -math.log10(pmin))
        s2 = 0.5  # 整合层内为常数（eQTL 存在是纳入判据）
        s3 = v['pip']
        sym_key = v['symbol'] or g
        s4 = norm(wd[sym_key], wmin, wmax) if sym_key in wd else 0.0
        score = 0.35 * s1 + 0.25 * s2 + 0.30 * s3 + 0.10 * s4
        score_rows.append({'gene': g, 'symbol': v['symbol'], 'score': round(score, 4),
                           'S1_DE': round(s1, 3), 'S2_eqtl': s2, 'S3_PIP': round(s3, 3),
                           'S4_network': round(s4, 3), 'n_queues': v['n_queues'],
                           'beta': v['beta'], 'p_DL': v['p_DL'], 'PIP': v['pip'],
                           'lead_SNP': v['snp'],
                           'STRING_degree': deg.get(sym_key, 0),
                           'STRING_weighted_degree': round(wd.get(sym_key, 0.0), 3),
                           'QTL_overlap': 'PENDING(QTLdb)'})
    score_rows.sort(key=lambda x: -x['score'])
    write_tsv('priority_scoring.tsv', list(score_rows[0].keys()), score_rows)

    def validation(r):
        """验证建议按可执行性分层（评审 M6）：无官方符号者不给 qPCR 引物建议。"""
        has_sym = bool(r['symbol'])
        if r['S3_PIP'] >= 0.9 and has_sym:
            return 'genotype lead SNP in mastitis-resistance cohorts; expression qPCR in mammary tissue'
        if r['S3_PIP'] >= 0.9 and not has_sym:
            return 'genotype lead SNP in resistance cohorts; transcript-model assay after annotation (Ensembl ID)'
        if has_sym:
            return 'expression qPCR panel (innate-immunity module); lead-SNP genotyping'
        return 'lead-SNP genotyping; annotation-dependent follow-up (no official symbol)'

    top = score_rows[:30]
    cols = ['rank', 'gene', 'symbol', 'score', 'S1_DE', 'S3_PIP', 'S4_network', 'n_queues',
            'beta', 'p_DL', 'lead_SNP', 'STRING_degree', 'STRING_weighted_degree',
            'validation_suggestion']
    rows = []
    for i, r in enumerate(top, 1):
        row = {'rank': i}
        for k in cols[1:-1]:
            row[k] = r[k]
        row['validation_suggestion'] = validation(r)
        rows.append(row)
    write_tsv('top_candidates.tsv', cols, rows)
    print('scored genes:', len(score_rows))
    for r in top[:10]:
        print(f"  {r['symbol'] or r['gene']:>12} score={r['score']:.3f} "
              f"(S1 {r['S1_DE']:.2f} S3 {r['S3_PIP']:.2f} S4 {r['S4_network']:.2f}) "
              f"pDL={r['p_DL']:.1e} deg={r['STRING_degree']}")
    write_json('priority_scoring_meta.json', {
        'n_scored': len(score_rows),
        'weights': {'S1_DE': 0.35, 'S2_eqtl': 0.25, 'S3_PIP': 0.30, 'S4_network': 0.10},
        'note': ('QTL layer pending (QTLdb maintenance); S4 uses the normalized STRING subnetwork '
                 'weighted degree; S2 is an inclusion criterion (constant within the integrated set); '
                 'STRING_degree is the unweighted subnetwork degree, STRING_weighted_degree the raw score sum')})


if __name__ == '__main__':
    main()
