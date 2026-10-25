#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 3a: 基因级合并 meta（芯片锚点基因 + RNA-seq 两部并入）。

合并口径：
  - 芯片轨：gene-anchored 1,886 基因（每基因取最显著探针的 β/p）
  - RNA-seq 轨：GSE131607（milk SC）+ GSE75379（udder 24h）welch t 结果（Ensembl 基因 ID 直配）
  - 跨队列合成同 04：Stouffer(√n) + DL 随机效应；RNA-seq 队列 n 小（6/6）权重自然低
  - 基因 ID 归一：Ensembl ENSBTAG 前缀去版本号；芯片锚点经 entrez_id 唯一化后无法直连
    Ensembl → 分两层报告（A: 全队列合成含 RNA-seq；B: 芯片轨基因级）
输出 work/output/de_meta_gene_combined.tsv
"""
import csv, gzip, json, os, re, sys
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, 'output')


def dl_random_effect(betas, ses):
    betas = np.asarray(betas, float); ses = np.asarray(ses, float)
    ok = ~np.isnan(betas) & ~np.isnan(ses)
    betas, ses = betas[ok], ses[ok]
    k = len(betas)
    if k == 0:
        return (np.nan,) * 4 + (0,)
    w = 1 / ses ** 2
    b_fe = (w * betas).sum() / w.sum()
    Q = (w * (betas - b_fe) ** 2).sum()
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max((Q - (k - 1)) / c, 0.0) if c > 0 else 0.0
    w_re = 1 / (ses ** 2 + tau2)
    b = (w_re * betas).sum() / w_re.sum()
    se = np.sqrt(1 / w_re.sum())
    z = b / se if se > 0 else np.nan
    return b, se, 2 * stats.norm.sf(abs(z)), Q, k


def se_from_bp(beta, p):
    p = min(max(p, 1e-300), 1 - 1e-16)
    z = stats.norm.isf(p / 2)
    return abs(beta) / z if z > 0 else np.nan


def main():
    # 1) 芯片轨基因级
    chip = {}
    for r in csv.DictReader(open(os.path.join(OUT, 'de_meta_gene_anchored.tsv'), encoding='utf-8'), delimiter='\t'):
        g = r['gene_symbol']
        chip.setdefault(g, []).append({'beta': float(r['beta_meta']), 'p': float(r['p_DL']),
                                        'n_q': int(r['n_queues'])})
    # 2) RNA-seq 轨 welch t（重算自 geo_parsed 的 z 矩阵 + 设计表分组）
    # GSE131607: 列名 H_/M_ 前缀
    def load_expr(gse):
        rows = list(csv.reader(gzip.open(os.path.join(HERE, 'data', 'geo_parsed', f'{gse}_expr.tsv.gz'), 'rt', encoding='utf-8'), delimiter='\t'))
        return rows[0][1:], [r[0] for r in rows[1:]], np.array([[float(x) for x in r[1:]] for r in rows[1:]])
    def welch_cols(mat, cols, case_pred, ctl_pred):
        ci = [i for i, c in enumerate(cols) if case_pred(c)]
        ki = [i for i, c in enumerate(cols) if ctl_pred(c)]
        x1, x0 = mat[:, ci], mat[:, ki]
        se = np.sqrt(x1.var(1, ddof=1) / len(ci) + x0.var(1, ddof=1) / len(ki))
        se[se == 0] = np.nan
        t = (x1.mean(1) - x0.mean(1)) / se
        p = 2 * stats.t.sf(np.abs(t), len(ci) + len(ki) - 2)
        return x1.mean(1) - x0.mean(1), p, len(ci), len(ki)

    cols, genes, mat = load_expr('GSE131607')
    b, p, n1, n0 = welch_cols(mat, cols, lambda c: c.startswith('M_'), lambda c: c.startswith('H_'))
    rnaseq = {'GSE131607': dict(zip(genes, zip(b, p, [n1 + n0] * len(genes))))}
    cols, genes, mat = load_expr('GSE75379')
    raw = gzip.open(os.path.join(HERE, 'data', 'geo', 'GSE75379_series_matrix.txt.gz'), 'rt', encoding='utf-8', errors='replace').read()
    titles = [v.strip('"') for v in re.search(r'^!Sample_title\s+(.*)$', raw, re.M).group(1).split('\t')]
    gsm = [v.strip('"') for v in re.search(r'^!Sample_geo_accession\s+(.*)$', raw, re.M).group(1).split('\t')]
    def char(field):
        for cm in re.finditer(r'^!Sample_characteristics_ch1\s+(.*)$', raw, re.M):
            cs = [v.strip().strip('"') for v in cm.group(1).split('\t')]
            if any(field in c.lower() for c in cs if ':' in c):
                return dict(zip(gsm, [c.split(':', 1)[1].strip() for c in cs]))
        return {}
    tissue, infect, timeh = char('tissue'), char('infection'), char('time (h)')
    title2gsm = dict(zip(titles, gsm))
    case_i, ctl_i = [], []
    for j, c in enumerate(cols):
        gsm0 = title2gsm.get(c.replace('_', '-'))
        if not gsm0:
            continue
        ti, inf, th = tissue.get(gsm0, ''), infect.get(gsm0, ''), timeh.get(gsm0, '')
        if ti.lower().startswith('udder') and th == '24':
            if inf.lower() == 'e.coli':
                case_i.append(j)
            elif 'kontrol' in inf.lower() or inf.lower() == 'control':
                ctl_i.append(j)
    b, p, n1, n0 = welch_cols(mat, cols, lambda j: False, lambda j: False)  # placeholder
    # 直接用索引计算
    x1, x0 = mat[:, case_i], mat[:, ctl_i]
    se = np.sqrt(x1.var(1, ddof=1) / len(case_i) + x0.var(1, ddof=1) / len(ctl_i))
    se[se == 0] = np.nan
    t = (x1.mean(1) - x0.mean(1)) / se
    p = 2 * stats.t.sf(np.abs(t), len(case_i) + len(ctl_i) - 2)
    bb = x1.mean(1) - x0.mean(1)
    rnaseq['GSE75379'] = dict(zip(genes, zip(bb, p, [len(case_i) + len(ctl_i)] * len(genes))))
    print('RNA-seq queues loaded:', {k: len(v) for k, v in rnaseq.items()})

    # 3) 合成：芯片基因经 entrez→Ensembl 桥接归一到 Ensembl ID 后与 RNA-seq 并入
    bridge = json.load(open(os.path.join(HERE, 'data', 'platform', 'entrez2ensembl_bridge.json')))
    anchor = {}
    for r in csv.DictReader(open(os.path.join(OUT, 'de_meta_gene_anchored.tsv'), encoding='utf-8'), delimiter='	'):
        e = r['entrez_id'].strip()
        ens = bridge.get(e)
        if not ens:
            continue
        anchor.setdefault(ens, []).append({'beta': float(r['beta_meta']), 'p': float(r['p_DL']),
                                            'n_q': int(r['n_queues']), 'sym': r['gene_symbol']})
    all_genes = set(anchor) | set(rnaseq['GSE131607']) | set(rnaseq['GSE75379'])
    rows_out = []
    sym_of = {}
    for g in all_genes:
        betas, ses, ns, srcs = [], [], [], []
        if g in anchor:
            best = sorted(anchor[g], key=lambda x: x['p'])[0]
            s = se_from_bp(best['beta'], best['p'])
            betas.append(best['beta']); ses.append(s); ns.append(0); srcs.append('chip')
            sym_of[g] = best['sym']
        for q in ('GSE131607', 'GSE75379'):
            if g in rnaseq[q]:
                bb, pp, nn = rnaseq[q][g]
                s = se_from_bp(float(bb), float(pp))
                betas.append(float(bb)); ses.append(s); ns.append(nn); srcs.append(q)
        if len(betas) < 2:
            continue  # 合成需 >=2 队列（单队列留在分轨结果里）
        b_re, se_re, p_re, Q, k = dl_random_effect(betas, ses)
        z_st = sum((be / s) * np.sqrt(n) for be, s, n in zip(betas, ses, ns) if s and s > 0 and n) , 0
        rows_out.append({'gene': g, 'symbol': sym_of.get(g, ''), 'n_queues': k,
                         'beta_meta': b_re, 'p_DL': p_re, 'Q': Q, 'sources': ','.join(srcs)})
    rows_out.sort(key=lambda r: r['p_DL'] if r['p_DL'] else 1)
    with open(os.path.join(OUT, 'de_meta_gene_combined.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gene', 'symbol', 'n_queues', 'beta_meta', 'p_DL', 'Q', 'sources'], delimiter='\t')
        w.writeheader()
        for r in rows_out:
            w.writerow({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items()})
    n2 = sum(1 for r in rows_out if r['n_queues'] >= 2)
    print(f'combined gene meta: {len(rows_out)} genes with >=2 queues; '
          f'p_DL<0.05: {sum(1 for r in rows_out if r["p_DL"] < 0.05)}')
    for r in rows_out[:12]:
        print(f"  {r['gene']:>16} q={r['n_queues']} β={r['beta_meta']:+.2f} p={r['p_DL']:.1e} [{r['sources']}]")


if __name__ == '__main__':
    main()
