#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 7c: 精确方差（exact Welch-t）敏感性分析（S10 的核心）。

评审 M7 指出主 meta 的 SE 由 p 的正态分位数反推（SE=|β|/z），而精确的每队列方差
应为 SE=|β|/|t|（t 为 Welch 统计量，df=n1+n0-2）。本脚本在全链上重算：

  1. 探针层：五队列（GPL2112）精确 SE → DL → BH → Tier A；与已发布口径对比
  2. 基因层：精确 p 重选锚点（每基因取最小 p 探针）→ 与两条 RNA-seq 轨（精确 SE）
     按 08 号完全相同的规则合成 → BH
  3. 评估既有 3,395 / 286 / Top-30 在精确口径下的存活：
     p<0.05、q<0.05、方向一致性、Top-30 排序重叠
  4. 报告 Tier A 的 I²（精确 SE 下）

输出 work/output/exact_gene_meta.tsv + exact_variance_sensitivity.json
运行：python work/code/20_exact_variance_sensitivity.py
"""
import csv, gzip, json, math, os, re, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.normpath(os.path.join(HERE, 'output'))
Q5 = ['GSE24560', 'GSE25319', 'GSE24217', 'GSE93082', 'GSE15025']


def outpath(name):
    base = os.path.normpath(OUT)
    p = os.path.normpath(base + os.sep + name)
    if not p.startswith(base + os.sep):
        raise SystemExit('output path escaped OUT dir: ' + name)
    return p


def dump(obj, name):
    Path(outpath(name)).write_text(json.dumps(obj, indent=1, ensure_ascii=False), encoding='utf-8')


def rd(name):
    return list(csv.DictReader(open(outpath(name), encoding='utf-8'), delimiter='\t'))


def rd_path(path):
    return list(csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'))


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    q = [1.0] * n
    run = 1.0
    for rank, i in enumerate(reversed(order), 1):
        run = min(run, ps[i] * n / max(n - rank + 1, 1))
        q[i] = run
    return q


def dl(betas, ses):
    b = np.asarray(betas, float)
    s = np.asarray(ses, float)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    b, s = b[ok], s[ok]
    k = len(b)
    if k == 0:
        return dict(k=0, beta=np.nan, se=np.nan, p=np.nan, Q=np.nan, I2=np.nan, tau2=0.0)
    w = 1 / s ** 2
    b_fe = (w * b).sum() / w.sum()
    Q = float((w * (b - b_fe) ** 2).sum())
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max((Q - (k - 1)) / c, 0.0) if c > 0 else 0.0
    w_re = 1 / (s ** 2 + tau2)
    beta = float((w_re * b).sum() / w_re.sum())
    se = float(np.sqrt(1 / w_re.sum()))
    p = float(2 * stats.norm.sf(abs(beta / se))) if se > 0 else np.nan
    I2 = max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0.0
    return dict(k=k, beta=beta, se=se, p=p, Q=Q, I2=float(I2), tau2=float(tau2))


def hksj(betas, ses, tau2):
    """Hartung-Knapp-Sidik-Jonkman：以 DL 的 tau2 为随机效应方差，t_{k-1} 推断。"""
    b = np.asarray(betas, float)
    s = np.asarray(ses, float)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    b, s = b[ok], s[ok]
    k = len(b)
    if k < 2:
        return np.nan
    w = 1 / (s ** 2 + tau2)
    b_re = (w * b).sum() / w.sum()
    var_hk = (1 / (k - 1)) * float((w * (b - b_re) ** 2).sum() / w.sum())
    if var_hk <= 0:
        return np.nan
    return float(2 * stats.t.sf(abs(b_re) / math.sqrt(var_hk), k - 1))


def load_expr(gse):
    rows = list(csv.reader(gzip.open(os.path.join(HERE, 'data', 'geo_parsed', f'{gse}_expr.tsv.gz'), 'rt', encoding='utf-8'), delimiter='\t'))
    return rows[0][1:], [r[0] for r in rows[1:]], np.array([[float(x) for x in r[1:]] for r in rows[1:]])


def welch_exact(x1, x0):
    n1, n0 = x1.shape[1], x0.shape[1]
    m1, m0 = x1.mean(1), x0.mean(1)
    se = np.sqrt(x1.var(1, ddof=1) / n1 + x0.var(1, ddof=1) / n0)
    se[se == 0] = np.nan
    t = (m1 - m0) / se
    df = n1 + n0 - 2
    p = 2 * stats.t.sf(np.abs(t), df)
    return m1 - m0, se, p, n1 + n0


def main():
    # ---------- 1. 探针层精确 SE ----------
    perq = defaultdict(dict)
    for r in rd('de_meta_chip_per_queue.tsv'):
        if r['p'] in ('NA', '') or r['t'] in ('NA', ''):
            continue
        perq[r['probe']][r['gse']] = (float(r['beta']), float(r['t']), float(r['p']))
    probe_rows = []
    for pb, d in perq.items():
        items = [(g, d[g]) for g in Q5 if g in d]
        if len(items) < 3:
            continue
        b = [x[1][0] for x in items]
        se = [abs(x[1][0] / x[1][1]) if x[1][1] != 0 else np.nan for x in items]
        res = dl(b, se)
        pos = sum(1 for x in b if x > 0)
        probe_rows.append({'probe': pb, 'beta': res['beta'], 'se': res['se'], 'p': res['p'],
                           'I2': res['I2'], 'p_hksj': hksj(b, se, res['tau2']),
                           'consist': max(pos, len(b) - pos) / len(b), 'k': res['k']})
    pq = bh([min(r['p'], 1.0) if math.isfinite(r['p']) else 1.0 for r in probe_rows])
    for r, q in zip(probe_rows, pq):
        r['q'] = q
        r['tierA'] = bool(q < 0.05 and r['consist'] == 1.0)
    pub_rows = rd('de_meta_chip_final.tsv')
    pub = {r['probe']: float(r['p_DL']) for r in pub_rows}
    pubA = {r['probe'] for r in pub_rows if r['tier'] == 'A'}
    exA = {r['probe'] for r in probe_rows if r['tierA']}
    common = [r for r in probe_rows if r['probe'] in pub]
    rho = stats.spearmanr([-math.log10(max(pub[r['probe']], 1e-300)) for r in common],
                          [-math.log10(max(r['p'], 1e-300)) for r in common]).statistic
    i2A = [r['I2'] for r in probe_rows if r['probe'] in pubA]
    hkA = sum(1 for r in probe_rows if r['probe'] in pubA and math.isfinite(r['p_hksj']) and r['p_hksj'] < 0.05)
    print(f'exact probe meta: n={len(probe_rows)} q<0.05={sum(1 for r in probe_rows if r["q"]<0.05)} tierA={len(exA)}')
    print(f'  published Tier A retained: {len(pubA & exA)}/{len(pubA)}; rho={rho:.4f}; I2 median={np.median(i2A):.1f}')

    # ---------- 2. 基因层精确 ----------
    plat = {}
    for r in rd_path(os.path.join(HERE, 'data', 'platform', 'GPL2112_probe2gene_final.tsv')):
        if r.get('gene_symbol'):
            plat.setdefault(r['probe'], r['gene_symbol'])
    best = {}
    for r in probe_rows:
        sym = plat.get(r['probe'])
        if not sym:
            continue
        cur = best.get(sym)
        if cur is None or (math.isfinite(r['p']) and r['p'] < cur['p']):
            best[sym] = {'beta': r['beta'], 'se': r['se'], 'p': r['p']}
    bridge = json.load(open(os.path.join(HERE, 'data', 'platform', 'entrez2ensembl_bridge.json'), encoding='utf-8'))
    sym2ens = {}
    for r in rd('de_meta_gene_anchored.tsv'):
        sym, e = r['gene_symbol'], r['entrez_id'].strip()
        if sym and e in bridge:
            sym2ens.setdefault(sym, bridge[e])
    # RNA-seq 精确
    rn = {}
    cols, genes, mat = load_expr('GSE131607')
    ci = [i for i, c in enumerate(cols) if c.startswith('M_')]
    ki = [i for i, c in enumerate(cols) if c.startswith('H_')]
    b, se, p, nn = welch_exact(mat[:, ci], mat[:, ki])
    rn['GSE131607'] = {g: (float(b[i]), float(se[i]), float(p[i]), int(nn)) for i, g in enumerate(genes)}
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
    t2g = dict(zip(titles, gsm))
    case_i, ctl_i = [], []
    for j, c in enumerate(cols):
        g0 = t2g.get(c.replace('_', '-'))
        if not g0:
            continue
        if tissue.get(g0, '').lower().startswith('udder') and timeh.get(g0) == '24':
            if infect.get(g0, '').lower() == 'e.coli':
                case_i.append(j)
            elif 'kontrol' in infect.get(g0, '').lower() or infect.get(g0, '').lower() == 'control':
                ctl_i.append(j)
    b, se, p, nn = welch_exact(mat[:, case_i], mat[:, ctl_i])
    rn['GSE75379'] = {g: (float(b[i]), float(se[i]), float(p[i]), int(nn)) for i, g in enumerate(genes)}

    # 合成（与 08 相同：>=2 个来源；DL）
    ens_betas = defaultdict(list)
    ens_ses = defaultdict(list)
    ens_syms = {}
    for sym, d in best.items():
        ens = sym2ens.get(sym)
        if not ens:
            continue
        ens_betas[ens].append(d['beta'])
        ens_ses[ens].append(d['se'])
        ens_syms[ens] = sym
    for q, table in rn.items():
        for g, (bb, ss, pp, nq) in table.items():
            ens_betas[g].append(bb)
            ens_ses[g].append(ss)
    gene_rows = []
    for g, bs in ens_betas.items():
        ss = ens_ses[g]
        if len(bs) < 2:
            continue
        res = dl(bs, ss)
        gene_rows.append({'gene': g, 'symbol': ens_syms.get(g, ''),
                          'n_sources': len(bs), 'beta': res['beta'], 'p': res['p'], 'q': None})
    qg = bh([min(r['p'], 1.0) if math.isfinite(r['p']) else 1.0 for r in gene_rows])
    for r, q in zip(gene_rows, qg):
        r['q'] = q
    exact = {r['gene']: r for r in gene_rows}
    Path(outpath('exact_gene_meta.tsv')).write_text(
        '\t'.join(['gene', 'symbol', 'n_sources', 'beta', 'p', 'q']) + '\n'
        + '\n'.join(f"{r['gene']}\t{r['symbol']}\t{r['n_sources']}\t{r['beta']:.6g}\t{r['p']:.6g}\t{r['q']:.6g}" for r in gene_rows)
        + '\n', encoding='utf-8')
    print(f'exact gene meta: n={len(gene_rows)} p<0.05={sum(1 for r in gene_rows if r["p"]<0.05)} q<0.05={sum(1 for r in gene_rows if r["q"]<0.05)}')

    # ---------- 3. 既有层在精确口径下的存活 ----------
    integ = rd('de_eqtl_pip_integrated.tsv')
    hi = [r for r in integ if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05]
    def survive(rows):
        out = {'n': len(rows), 'with_exact_p': 0, 'p_lt_0.05': 0, 'q_lt_0.05': 0, 'sign_same': 0}
        for r in rows:
            e = exact.get(r['gene'])
            if not e:
                continue
            out['with_exact_p'] += 1
            if math.isfinite(e['p']) and e['p'] < 0.05:
                out['p_lt_0.05'] += 1
            if math.isfinite(e['q']) and e['q'] < 0.05:
                out['q_lt_0.05'] += 1
            if float(r['beta_meta']) * e['beta'] > 0:
                out['sign_same'] += 1
        return out
    hi_exact_p = sum(1 for r in hi if r['gene'] in exact and exact[r['gene']]['p'] < 0.05)
    # 精确口径下的 286 类比层：PIP>=0.5 且有精确 p<0.05 / q<0.05
    hi_like = sum(1 for r in integ if float(r['max_PIP']) >= 0.5 and r['gene'] in exact and exact[r['gene']]['p'] < 0.05)
    hi_like_q = sum(1 for r in integ if float(r['max_PIP']) >= 0.5 and r['gene'] in exact and exact[r['gene']]['q'] < 0.05)

    # ---------- 4. Top-30 在精确口径下的排序 ----------
    sc = rd('priority_scoring.tsv')
    S3 = {r['gene']: float(r['S3_PIP']) for r in sc}
    S4 = {r['gene']: float(r['S4_network']) for r in sc}
    base_top30 = {r['gene'] for r in sc[:30]}
    lp = [(-math.log10(max(exact[g]['p'], 1e-300)), g) for g in exact if g in S3 and math.isfinite(exact[g]['p'])]
    if lp:
        lo = min(x[0] for x in lp)
        hii = max(x[0] for x in lp)
        rng = (hii - lo) or 1.0
        scores = []
        for val, g in lp:
            s1 = (val - lo) / rng
            scores.append((0.35 * s1 + 0.25 * 0.5 + 0.30 * S3[g] + 0.10 * S4.get(g, 0.0), g))
        scores.sort(reverse=True)
        exact_top30 = {g for _, g in scores[:30]}
        ov = len(exact_top30 & base_top30)
        exact_top10 = {g for _, g in scores[:10]}
        ov10 = len(exact_top10 & {r['gene'] for r in sc[:10]})
    else:
        ov = ov10 = None

    res = {
        'probe_level': {
            'probes_exact': len(probe_rows), 'q_lt_0.05_exact': sum(1 for r in probe_rows if r['q'] < 0.05),
            'tierA_exact': len(exA), 'published_tierA': len(pubA),
            'published_tierA_retained_exact': len(pubA & exA),
            'published_tierA_retained_pct': round(100 * len(pubA & exA) / len(pubA), 1),
            'spearman_rho_logp': round(float(rho), 4),
            'tierA_I2_median_exact_pct': round(float(np.median(i2A)), 1) if i2A else None,
            'tierA_hksj_p05': hkA,
        },
        'gene_level': {
            'genes_exact': len(gene_rows), 'p_lt_0.05': sum(1 for r in gene_rows if r['p'] < 0.05),
            'q_lt_0.05': sum(1 for r in gene_rows if r['q'] < 0.05),
            'integrated_3395': survive(integ),
            'high_conf_286': survive(hi),
            'high_conf_286_exact_p_lt_0.05': hi_exact_p,
            'high_conf_like_exact_p_lt_0.05': hi_like,
            'high_conf_like_exact_q_lt_0.05': hi_like_q,
            'top30_overlap_exact': ov, 'top10_overlap_exact': ov10,
        },
        'note': 'exact SE = |beta/t| from the per-queue Welch statistic; gene-level rebuild replicates the stage-3 rules (>=2 sources, DL)',
    }
    dump(res, 'exact_variance_sensitivity.json')
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == '__main__':
    main()
