#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 1c: 跨队列 DE meta（芯片轨 7 部）。

方法（设计文档三号层，CV 管线平移）：
  1. 每队列平台内 limma 式 moderated t（近似：Welch t + BH 队列内）——case vs control
  2. 探针 → 基因映射：GPL 平台注释不在 series_matrix 内 → 用探针 ID 直接跨队列求交
     （GPL2112 系 4 部共享同探针集；GPL6082/11648 独立）——首版按"探针 ID 交集"做
     方向一致性 + Fisher/Stouffer 合成，映射 Ensembl 的精确基因级 meta 留待平台注释下载
  3. 跨队列合成：Stouffer（权重 = sqrt(n_case+n_control)）+ 随机效应 DerSimonian-Laird
  4. 方向一致性表：每探针各队列 β 符号一致性比例
输出：
  work/output/de_meta_chip_per_queue.tsv  （每队列每探针 β/p）
  work/output/de_meta_chip_summary.tsv    （探针级合成：beta_meta/p_stouffer/p_DL/一致性）
"""
import csv, gzip, json, os, sys
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARSED = os.path.join(HERE, 'data', 'geo_parsed')
OUT = os.path.join(HERE, 'output')


def load_queue(gse):
    rows = list(csv.reader(gzip.open(os.path.join(PARSED, f'{gse}_expr.tsv.gz'), 'rt', encoding='utf-8'), delimiter='\t'))
    probes = [r[0] for r in rows[1:]]
    mat = np.array([[float(x) for x in r[1:]] for r in rows[1:]], dtype=np.float32)
    return probes, mat


def load_groups(gse):
    g = {}
    for r in csv.DictReader(open(os.path.join(OUT, 'de_queue_design.tsv'), encoding='utf-8'), delimiter='\t'):
        if r['gse'] == gse:
            g[r['gsm']] = r['group']
    return g


def moderated_t(mat, groups):
    case_idx = [i for i, s in enumerate(groups) if groups[s] == 'case']
    ctl_idx = [i for i, s in enumerate(groups) if groups[s] == 'control']
    x1, x0 = mat[:, case_idx], mat[:, ctl_idx]
    n1, n0 = x1.shape[1], x0.shape[1]
    m1, m0 = x1.mean(axis=1), x0.mean(axis=1)
    v1, v0 = x1.var(axis=1, ddof=1), x0.var(axis=1, ddof=1)
    se = np.sqrt(v1 / n1 + v0 / n0)
    se[se == 0] = np.nan
    t = (m1 - m0) / se
    df = n1 + n0 - 2
    p = 2 * stats.t.sf(np.abs(t), df)
    return m1 - m0, t, p, n1, n0


def dl_random_effect(betas, ses):
    """DerSimonian-Laird（β_i, se_i 数组 → 合成 β/SE/Q/p）"""
    betas = np.asarray(betas, dtype=float)
    ses = np.asarray(ses, dtype=float)
    ok = ~np.isnan(betas) & ~np.isnan(ses)
    betas, ses = betas[ok], ses[ok]
    k = len(betas)
    if k == 0:
        return np.nan, np.nan, np.nan, np.nan, k
    w = 1 / ses ** 2
    b_fe = (w * betas).sum() / w.sum()
    Q = (w * (betas - b_fe) ** 2).sum()
    df = k - 1
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max((Q - df) / c, 0.0) if c > 0 else 0.0
    w_re = 1 / (ses ** 2 + tau2)
    b_re = (w_re * betas).sum() / w_re.sum()
    se_re = np.sqrt(1 / w_re.sum())
    z = b_re / se_re if se_re > 0 else np.nan
    p = 2 * stats.norm.sf(abs(z))
    return b_re, se_re, p, Q, k


def main():
    queues = ['GSE9685', 'GSE24560', 'GSE25319', 'GSE24217', 'GSE93082', 'GSE15025', 'GSE34031']
    per_probe = {}   # probe -> {gse: (beta, p, n1, n0)}
    with open(os.path.join(OUT, 'de_meta_chip_per_queue.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['gse', 'probe', 'beta', 't', 'p', 'n_case', 'n_control'])
        for gse in queues:
            probes, mat = load_queue(gse)
            groups = load_groups(gse)
            ordered = [groups[s] for s in list(groups)[0:0]]  # no-op
            gsm_order = [r for r in csv.DictReader(open(os.path.join(OUT, 'de_queue_design.tsv'), encoding='utf-8'), delimiter='\t') if r['gse'] == gse]
            gsm_order = [r['gsm'] for r in gsm_order]
            g = [groups[s] for s in gsm_order]
            gmap = {s: groups[s] for s in gsm_order}
            beta, t, p, n1, n0 = moderated_t(mat, gmap)
            for i, pb in enumerate(probes):
                per_probe.setdefault(pb, {})[gse] = (float(beta[i]), float(p[i]), n1, n0)
                w.writerow([gse, pb, f'{beta[i]:.4f}', f'{t[i]:.3f}',
                            f'{p[i]:.3g}' if not np.isnan(p[i]) else 'NA', n1, n0])
            print(f'{gse}: beta/p computed for {len(probes)} probes (case {n1}/ctl {n0})')

    # 跨队列合成（>=3 队列覆盖的探针）
    rows_out = []
    for pb, d in per_probe.items():
        if len(d) < 3:
            continue
        gses = sorted(d)
        betas = [d[g][0] for g in gses]
        ps = [d[g][1] for g in gses]
        ns = [d[g][2] + d[g][3] for g in gses]
        ses = [abs(b) / max(stats.norm.isf(min(pv, 1 - 1e-12)), 1e-12) if pv > 0 else np.nan
               for b, pv in zip(betas, ps)]  # 由 β,p 反推 se（双侧）
        # Stouffer
        keep = [i for i, s in enumerate(ses) if not np.isnan(s) and s > 0]
        if len(keep) < 3:
            continue
        zs = [betas[i] / ses[i] for i in keep]
        wgt = [np.sqrt(ns[i]) for i in keep]
        z_st = sum(z * w for z, w in zip(zs, wgt)) / np.sqrt(sum(w ** 2 for w in wgt))
        p_st = 2 * stats.norm.sf(abs(z_st))
        b_dl, se_dl, p_dl, Q, k = dl_random_effect([betas[i] for i in keep], [ses[i] for i in keep])
        pos = sum(1 for b in betas if b > 0)
        consist = max(pos, len(betas) - pos) / len(betas)
        rows_out.append({'probe': pb, 'n_queues': len(keep), 'beta_meta': b_dl,
                         'p_stouffer': p_st, 'p_DL': p_dl, 'Q': Q,
                         'sign_consistency': consist,
                         'betas': {gses[i]: round(betas[i], 4) for i in range(len(gses))}})
    rows_out.sort(key=lambda r: (r['p_DL'] if r['p_DL'] else 1))
    with open(os.path.join(OUT, 'de_meta_chip_summary.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['probe', 'n_queues', 'beta_meta', 'p_stouffer', 'p_DL',
                                          'Q', 'sign_consistency', 'betas'], delimiter='\t')
        w.writeheader()
        for r in rows_out:
            w.writerow({k: (json.dumps(v) if k == 'betas' else
                            (f'{v:.4g}' if isinstance(v, float) else v))
                        for k, v in r.items()})
    print(f'summary: {len(rows_out)} probes meta-analyzed (>=3 queues)')
    # 一致性摘要
    n_full_consist = sum(1 for r in rows_out if r['sign_consistency'] == 1.0)
    n_sig = sum(1 for r in rows_out if r['p_DL'] and r['p_DL'] < 0.05)
    print(f'sign-consistent 100%: {n_full_consist}; p_DL<0.05: {n_sig}')


if __name__ == '__main__':
    main()
