#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 7b: 评审响应补充分析（S10 的数据源）。

针对模拟同行评审 M1/M2/M4/M5/M7/M10 与次要项，逐项给出可复算数字：

  A 基因层 BH-FDR（combined / integrated / 286 层 / Top-30）与期望假阳性
  B 五个已知轴基因在 combined / eQTL-overlap / integrated 中的位置（回答"为何不进 Top-30"）
  C 竞文 150 标识的符号级重映射（Ensembl REST，bos_taurus）与重算交集（M5）
  D Danish 肝臂 vs 乳臂效应估计相关（M4 独立性量化）
  E SE 由 Welch t 直接导出重算（M7）+ I²/τ² + Hartung-Knapp-Sidik-Jonkman 敏感性
  F Top-30 权重网格稳定性（M2/次要 16）
  G 五个贡献队列的缺失率与 Tier A 覆盖（次要 11）
  H 探针层 Tier A 上/下调拆分（次要 7）
  I 队列贡献账（9 队列 → 定量合成 7 队列 / 257 样本）

输出 work/output/review_response_metrics.json + competitor_symbol_mapping.json
运行：在 VM2 课题根目录 python work/code/19_review_response.py
"""
import csv, json, math, os, ssl, sys, time, urllib.parse, urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.normpath(os.path.join(HERE, 'output'))
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}), urllib.request.HTTPSHandler(context=CTX))
OPENER.addheaders = [('User-Agent', 'VM2-review/1.0')]

CHIP_CONTRIB = ['GSE24560', 'GSE25319', 'GSE24217', 'GSE93082', 'GSE15025']


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


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    q = [1.0] * n
    run = 1.0
    for rank, i in enumerate(reversed(order), 1):
        run = min(run, ps[i] * n / max(n - rank + 1, 1))
        q[i] = run
    return q


def finite(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def dl(betas, ses):
    b = np.asarray(betas, float)
    s = np.asarray(ses, float)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    b, s = b[ok], s[ok]
    k = len(b)
    if k == 0:
        return dict(k=0, beta=np.nan, se=np.nan, p=np.nan, Q=np.nan, tau2=np.nan, I2=np.nan)
    w = 1 / s ** 2
    b_fe = (w * b).sum() / w.sum()
    Q = float((w * (b - b_fe) ** 2).sum())
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max((Q - (k - 1)) / c, 0.0) if c > 0 else 0.0
    w_re = 1 / (s ** 2 + tau2)
    beta = float((w_re * b).sum() / w_re.sum())
    se = float(np.sqrt(1 / w_re.sum()))
    p = float(2 * stats.norm.sf(abs(beta / se))) if se > 0 else np.nan
    I2 = (Q - (k - 1)) / Q * 100 if Q > 0 else 0.0
    return dict(k=k, beta=beta, se=se, p=p, Q=Q, tau2=float(tau2), I2=float(max(0.0, I2)))


def hksj(betas, ses):
    """Hartung-Knapp-Sidik-Jonkman（随机效应，t_{k-1}）。"""
    b = np.asarray(betas, float)
    s = np.asarray(ses, float)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    b, s = b[ok], s[ok]
    k = len(b)
    if k < 2:
        return np.nan
    d = dl(b, s)
    w = 1 / (s ** 2 + d['tau2'])
    b_re = (w * b).sum() / w.sum()
    var_hk = (1 / (k - 1)) * float((w * (b - b_re) ** 2).sum() / w.sum())
    if var_hk <= 0:
        return np.nan
    return float(2 * stats.t.sf(abs(b_re) / math.sqrt(var_hk), k - 1))


def ensembl_lookup_symbol(sym):
    url = ('https://rest.ensembl.org/lookup/symbol/bos_taurus/'
           + urllib.parse.quote(sym) + '?content-type=application/json')
    try:
        d = json.loads(OPENER.open(url, timeout=40).read())
        return d.get('id'), d.get('assembly_name')
    except Exception:
        return None, None


def ensembl_lookup_id(gid):
    url = ('https://rest.ensembl.org/lookup/id/' + urllib.parse.quote(gid)
           + '?content-type=application/json')
    try:
        d = json.loads(OPENER.open(url, timeout=40).read())
        return d.get('display_name'), d.get('assembly_name')
    except Exception:
        return None, None


def main():
    res = {}

    comb = rd('de_meta_gene_combined.tsv')
    integ = rd('de_eqtl_pip_integrated.tsv')
    sc = rd('priority_scoring.tsv')
    fin = rd('de_meta_chip_final.tsv')
    perq = rd('de_meta_chip_per_queue.tsv')

    # ---------- A. 基因层 BH ----------
    vals = [(r, float(r['p_DL'])) for r in comb if finite(r['p_DL'])]
    qs = bh([p for _, p in vals])
    n_p05 = sum(1 for _, p in vals if p < 0.05)
    n_q05 = sum(1 for q in qs if q < 0.05)
    ivals = [(r, float(r['p_DL'])) for r in integ if finite(r['p_DL'])]
    iqs = bh([p for _, p in ivals])
    qmap = {r['gene']: q for (r, _), q in zip(ivals, iqs)}
    hi = [r for r in integ if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05]
    hi_q = sum(1 for r in hi if qmap.get(r['gene'], 1.0) < 0.05)
    top30 = sc[:30]
    top30_q = sum(1 for r in top30 if qmap.get(r['gene'], 1.0) < 0.05)
    res['A_gene_level_fdr'] = {
        'combined_genes_with_defined_p': len(vals),
        'combined_p_lt_0.05': n_p05,
        'combined_q_lt_0.05': n_q05,
        'combined_expected_fp_at_p05': round(0.05 * len(vals), 1),
        'integrated_genes_with_defined_p': len(ivals),
        'integrated_p_lt_0.05': sum(1 for _, p in ivals if p < 0.05),
        'integrated_q_lt_0.05': sum(1 for q in iqs if q < 0.05),
        'integrated_expected_fp_at_p05': round(0.05 * len(ivals), 1),
        'high_conf_286': len(hi),
        'high_conf_with_gene_q_lt_0.05': hi_q,
        'high_conf_expected_fp_at_p05': round(0.05 * len(hi), 1),
        'top30_with_gene_q_lt_0.05': top30_q,
        'note': 'gene-level BH over genes with a defined combined (DL) p; the 286 tier is a PIP-based screening tier, not an FDR-controlled tier',
    }

    # ---------- B. 轴基因定位 ----------
    sym2row = {}
    for r in comb:
        if r.get('symbol'):
            sym2row.setdefault(r['symbol'], []).append(r)
    overlap_ens = set()
    op = Path(outpath('de_eqtl_overlap_genes.txt'))
    if op.is_file():
        overlap_ens = set(line.split('\t')[0].strip() for line in op.read_text(encoding='utf-8').splitlines()
                          if line.strip())
    integ_ens = {r['gene'] for r in integ}
    axis = {}
    for s in ('SAA3', 'CCL2', 'SOCS3', 'TIMP1', 'SERPINE1', 'RNASE1', 'FOS', 'NFKBIA', 'PRL', 'MET', 'CD3G'):
        rows = sym2row.get(s, [])
        row = rows[0] if rows else None
        ens = row['gene'] if row else None
        rank = None
        if row and finite(row['p_DL']):
            srt = sorted([r for r in comb if finite(r['p_DL'])], key=lambda r: float(r['p_DL']))
            rank = next((i + 1 for i, r in enumerate(srt) if r is row), None)
        axis[s] = {
            'ensembl': ens, 'beta': row['beta_meta'] if row else None,
            'p_DL': row['p_DL'] if row else None, 'rank_by_p': rank,
            'in_eqtl_overlap_6203': bool(ens and ens in overlap_ens),
            'in_integrated_3395': bool(ens and ens in integ_ens),
            'pip': next((r['max_PIP'] for r in integ if r['gene'] == ens), None),
        }
    res['B_axis_genes'] = axis

    # ---------- C. 竞文符号重映射 ----------
    comp = json.load(open(outpath('mini_qtl_gene_baseline.json'), encoding='utf-8'))
    ens_ids = [c for c in comp if c.startswith('ENSBTAG')]
    syms = [c for c in comp if not c.startswith('ENSBTAG')]
    mapping = {'source': 'Ensembl REST /lookup/symbol/bos_taurus (accessed 2026-10-25)',
               'n_identifiers': len(comp), 'n_ensembl': len(ens_ids), 'n_symbols': len(syms),
               'symbol_map': {}, 'unresolved': []}
    mapped_ens = set(ens_ids)
    for s in syms:
        gid, asm = ensembl_lookup_symbol(s)
        if gid:
            mapping['symbol_map'][s] = {'ensembl': gid, 'assembly': asm}
            mapped_ens.add(gid)
        else:
            mapping['unresolved'].append(s)
        time.sleep(0.05)
    shared = {}
    for gid in ens_ids[:0]:
        pass
    # shared gene symbol
    shared_ids = set(c['gene'] for c in [])
    hi_ens = {r['gene'] for r in hi}
    int_ens = {r['gene'] for r in integ}
    inter_hi = sorted(mapped_ens & hi_ens)
    inter_int = sorted(mapped_ens & int_ens)
    shared_symbols = {}
    for gid in sorted(set(inter_hi) | set(inter_int)):
        nm, asm = ensembl_lookup_id(gid)
        shared_symbols[gid] = {'symbol': nm, 'assembly': asm}
    mapping['intersection_high_conf_286'] = inter_hi
    mapping['intersection_integrated_3395'] = inter_int
    mapping['shared_symbols'] = shared_symbols
    mapping['n_symbol_resolved'] = len(mapping['symbol_map'])
    mapping['as_published_intersection'] = ['ENSBTAG00000035158']
    res['C_competitor_mapping'] = {
        'n_identifiers': len(comp), 'n_ensembl_given': len(ens_ids), 'n_symbols': len(syms),
        'n_symbols_resolved': len(mapping['symbol_map']), 'n_unresolved': len(mapping['unresolved']),
        'mapped_total_unique': len(mapped_ens),
        'intersection_286_as_published': 1,
        'intersection_286_symbol_resolved': len(inter_hi),
        'intersection_3395_symbol_resolved': len(inter_int),
        'shared_genes_in_286': shared_symbols,
    }
    dump(mapping, 'competitor_symbol_mapping.json')

    # ---------- D. Danish 双组织臂相关 ----------
    dq = defaultdict(dict)
    for r in perq:
        if r['p'] in ('NA', ''):
            continue
        dq[r['probe']][r['gse']] = (float(r['beta']), float(r['t']), float(r['p']))
    common = [p for p in dq if all(g in dq[p] for g in CHIP_CONTRIB)]
    x = np.array([dq[p]['GSE25319'][0] for p in common])
    y = np.array([dq[p]['GSE24217'][0] for p in common])
    r_arm = float(np.corrcoef(x, y)[0, 1])
    res['D_danish_arms'] = {
        'probes_all_five': len(common),
        'liver_25319_vs_udder_24217_beta_pearson_r': round(r_arm, 4),
        'note': 'same challenged cows, different tissues; partial correlation (not duplication)',
    }

    # ---------- E. SE 由 t 导出 + I²/τ² + HKSJ ----------
    pub = {}
    for r in fin:
        pb = r['probe']
        pub[pb] = dict(beta=float(r['beta_meta']), p=float(r['p_DL']), q=float(r['q_BH']), tier=r['tier'])
    probes_e = []
    for p in common:
        entries = [dq[p][g] for g in CHIP_CONTRIB if g in dq[p]]
        if len(entries) < 3:
            continue
        betas = [e[0] for e in entries]
        ts = [e[1] for e in entries]
        ses_t = [abs(b / t) if t not in (0.0,) else np.nan for b, t in zip(betas, ts)]
        ses_z = [abs(b) / max(stats.norm.isf(min(e[2], 1 - 1e-12)), 1e-12) if b != 0 else np.nan
                 for b, e in zip(betas, entries)]
        d_t = dl(betas, ses_t)
        d_z = dl(betas, ses_z)
        probes_e.append({'probe': p, 'beta': d_t['beta'], 'p_t': d_t['p'], 'p_z': d_z['p'],
                         'I2': d_z['I2'], 'tau2': d_z['tau2'], 'k': d_t['k']})
    q_t = bh([min(x['p_t'], 1.0) if np.isfinite(x['p_t']) else 1.0 for x in probes_e])
    q_z = bh([min(x['p_z'], 1.0) if np.isfinite(x['p_z']) else 1.0 for x in probes_e])
    # Tier A (published) evaluation in the t-SE recomputation
    pubA = [x for x in probes_e if pub.get(x['probe'], {}).get('tier') == 'A']
    both_sig = sum(1 for x, qt in zip(probes_e, q_t) if pub.get(x['probe'], {}).get('tier') == 'A' and qt < 0.05)
    sign_same = sum(1 for x in pubA if pub[x['probe']]['beta'] * x['beta'] > 0)
    rho_pt = stats.spearmanr([-math.log10(max(pub[x['probe']]['p'], 1e-300)) for x in pubA],
                             [-math.log10(max(x['p_t'], 1e-300)) for x in pubA]).statistic
    # HKSJ on published Tier A
    hk_sig = 0
    hk_p = {}
    for x in pubA:
        entries = [dq[x['probe']][g] for g in CHIP_CONTRIB if g in dq[x['probe']]]
        if len(entries) < 3:
            continue
        betas = [e[0] for e in entries]
        ses = [abs(b) / max(stats.norm.isf(min(e[2], 1 - 1e-12)), 1e-12) if b != 0 else np.nan
               for b, e in zip(betas, entries)]
        ph = hksj(betas, ses)
        hk_p[x['probe']] = ph
        if np.isfinite(ph) and ph < 0.05:
            hk_sig += 1
    i2A = [x['I2'] for x in pubA if np.isfinite(x['I2'])]
    res['E_se_from_t'] = {
        'probes_analyzed': len(probes_e),
        'tierA_n': len(pubA),
        'tierA_still_q05_with_t_se': both_sig,
        'tierA_sign_concordance_pct': round(100 * sign_same / len(pubA), 2),
        'spearman_rho_logp_published_vs_t': round(float(rho_pt), 4),
        'q05_total_published_convention': sum(1 for qq in q_z if qq < 0.05),
        'q05_total_t_convention': sum(1 for qq in q_t if qq < 0.05),
        'tierA_I2_median_pct': round(float(np.median(i2A)), 1) if i2A else None,
        'tierA_I2_q1_q3': [round(float(np.percentile(i2A, 25)), 1), round(float(np.percentile(i2A, 75)), 1)] if i2A else None,
        'tierA_hksj_p05': hk_sig,
        'note': 't-SE is the exact inversion from the per-queue Welch t; published convention uses the normal quantile of the two-sided p',
        'superseded_by': 'exact_variance_sensitivity.json (20_exact_variance_sensitivity.py, >=3-queue coverage; this block is a five-queue-complete subset only)',
    }

    # ---------- F. 权重网格 ----------
    s1 = np.array([float(r['S1_DE']) for r in sc])
    s3 = np.array([float(r['S3_PIP']) for r in sc])
    s4 = np.array([float(r['S4_network']) for r in sc])
    base = 0.35 * s1 + 0.30 * s3 + 0.10 * s4
    top30 = set(np.argsort(-base)[:30])
    top10 = set(np.argsort(-base)[:10])
    grid = []
    step = 0.05
    for w1i in range(1, 14):
        for w3i in range(1, 14):
            w1, w3 = round(w1i * step, 2), round(w3i * step, 2)
            w4 = round(0.75 - w1 - w3, 2)
            if w4 < 0.05:   # symmetric grid: all three active weights >= 0.05
                continue
            scj = w1 * s1 + w3 * s3 + w4 * s4
            order = np.argsort(-scj)
            grid.append({'w_DE': w1, 'w_PIP': w3, 'w_network': w4,
                         'top30_overlap': len(set(order[:30]) & top30),
                         'top10_overlap': len(set(order[:10]) & top10)})
    ov = [g['top30_overlap'] for g in grid]
    res['F_weight_grid'] = {
        'grid_points': len(grid), 'step': step,
        'top30_overlap_median': float(np.median(ov)),
        'top30_overlap_min': int(np.min(ov)), 'top30_overlap_max': int(np.max(ov)),
        'grid_points_with_overlap_ge_27': sum(1 for v in ov if v >= 27),
        'grid_points_with_overlap_30': sum(1 for v in ov if v == 30),
        'top10_overlap_min': int(min(g['top10_overlap'] for g in grid)),
        'note': 'weights sum to 0.75 (eQTL-presence term fixed at 0.25, constant across scored genes)',
    }
    Path(outpath('weight_grid.tsv')).write_text(
        '\t'.join(['w_DE', 'w_PIP', 'w_network', 'top30_overlap', 'top10_overlap']) + '\n'
        + '\n'.join('\t'.join(str(g[k]) for k in ('w_DE', 'w_PIP', 'w_network', 'top30_overlap', 'top10_overlap')) for g in grid)
        + '\n', encoding='utf-8')

    # ---------- G. 缺失率（五个贡献队列） ----------
    import gzip
    tierA = {r['probe'] for r in fin if r['tier'] == 'A'}
    miss = {}
    for gse in CHIP_CONTRIB:
        path = Path(os.path.join(HERE, 'data', 'geo', gse + '_series_matrix.txt.gz'))
        if not path.is_file():
            continue
        frac = {}
        started = False
        hdr = None
        with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('"ID_REF"'):
                    hdr = line.rstrip('\n').split('\t')
                    started = True
                    continue
                if started:
                    if not line.strip() or line.startswith('!') or line.startswith('series_matrix_table_end'):
                        break
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 3:
                        continue
                    pid = parts[0].strip('"')
                    vals = parts[1:]
                    n = len(vals)
                    ok = sum(1 for v in vals if v.strip('"').strip() != '')
                    frac[pid] = ok / n if n else 0.0
        if not frac:
            continue
        ta = [frac[p] for p in frac if p in tierA]
        miss[gse] = {
            'probes': len(frac),
            'median_nonmissing_pct': round(100 * float(np.median(list(frac.values()))), 1),
            'pct_probes_ge80_nonmissing': round(100 * sum(1 for v in frac.values() if v >= 0.8) / len(frac), 1),
            'tierA_probes_measured': len(ta),
            'tierA_pct_ge80_nonmissing': round(100 * sum(1 for v in ta if v >= 0.8) / len(ta), 1) if ta else None,
        }
    res['G_missingness'] = miss

    # ---------- H. Tier A 上/下调 ----------
    up = sum(1 for r in fin if r['tier'] == 'A' and float(r['beta_meta']) > 0)
    down = sum(1 for r in fin if r['tier'] == 'A' and float(r['beta_meta']) < 0)
    au = sum(1 for r in rd('de_meta_gene_anchored.tsv') if float(r['beta_meta']) > 0)
    ad = sum(1 for r in rd('de_meta_gene_anchored.tsv') if float(r['beta_meta']) < 0)
    res['H_tierA_direction'] = {'probe_tierA_up': up, 'probe_tierA_down': down,
                                'chip_anchor_genes_up': au, 'chip_anchor_genes_down': ad}

    # ---------- I. 队列贡献账 ----------
    meta = json.load(open(outpath('de_meta_queue_metadata.json'), encoding='utf-8'))['queue']
    total = sum(int(v.get('n_samples') or 0) for v in meta.values())
    res['I_cohort_accounting'] = {
        'cohorts_assembled': len(meta), 'samples_assembled': total,
        'chip_queues_in_probe_meta': CHIP_CONTRIB,
        'chip_samples_in_probe_meta': 88 + 54 + 49 + 12 + 30,
        'rnaseq_cohorts_in_gene_meta': ['GSE131607', 'GSE75379'],
        'rnaseq_samples_in_gene_meta': 12 + 12,
        'samples_in_quantitative_synthesis': 233 + 24,
        'excluded_not_combinable': {'GSE9685': 36, 'GSE34031': 16,
                                    'reason': 'deposited series matrices carry numeric row indices, not probe identifiers'},
        'excluded_by_design': {'GSE75379_liver_arm': 24, 'GSE93082_endometritis': 6},
        'samples_with_per_queue_contrasts': 285 + 24,
    }

    dump(res, 'review_response_metrics.json')
    print(json.dumps(res, indent=1, ensure_ascii=False)[:6000])


if __name__ == '__main__':
    main()
