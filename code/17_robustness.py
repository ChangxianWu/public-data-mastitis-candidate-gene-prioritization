#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 2b: 稳健性三件套（S8 补充材料的数据源）。

  A. 病原分层探针级 meta（S3 分组口径）：E. coli {25319,24217,15025} vs
     S. aureus 相关 {9685,24560,93082}；field 层 {34031,131607} 探针空间不相交
     （Illumina + RNA-seq），不可合并，仅记录。
  B. 去肝臂敏感性：主 meta 探针池去 GSE25319 后重合成（肝臂/乳臂同批 Danish
     攻毒动物的整臂移除，量化 LOQO 单队列下探的机理）。
  C. Top-30 权重扰动稳定性：S1/S3/S4 权重各 ±50% 抖动重归一（S2 层为常数不参与
     排序），1,000 次 + 三个单分量归零压力测试。

方法与 04 号主 meta 逐位一致（同模块导入：Welch t、se 反推、Stouffer sqrt(n)、
DL 随机效应、层内 BH），层内覆盖规则同为 >=3 队列。
输出 work/output/stratified_meta.json + weight_perturbation.json。
"""
import csv, importlib.util, json, os, sys
from pathlib import Path
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.normpath(os.path.join(HERE, 'output'))


def outpath(name):
    """OUT 目录下的静态文件名拼接：归一化并校验不越出 OUT（安全门要求）。"""
    base = os.path.normpath(OUT)
    p = os.path.normpath(base + os.sep + name)
    if not p.startswith(base + os.sep):
        raise SystemExit('output path escaped OUT dir: ' + name)
    return p


def dump_json(obj, name):
    Path(outpath(name)).write_text(json.dumps(obj, indent=1), encoding='utf-8')


spec = importlib.util.spec_from_file_location('de_meta_04', os.path.join(HERE, 'code', '04_de_meta.py'))
m04 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m04)

CHIP_QUEUES = ['GSE9685', 'GSE24560', 'GSE25319', 'GSE24217', 'GSE93082', 'GSE15025', 'GSE34031']
ECOLI = ['GSE25319', 'GSE24217', 'GSE15025']
SAUR = ['GSE9685', 'GSE24560', 'GSE93082']
NOLIVER = ['GSE9685', 'GSE24560', 'GSE24217', 'GSE93082', 'GSE15025']


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    q = [1.0] * n
    running = 1.0
    for rank, i in enumerate(reversed(order), 1):
        running = min(running, ps[i] * n / (n - rank + 1))
        q[i] = running
    return q


def compute_per_queue():
    """与 04.main 逐位相同的每队列 Welch（同一 design 行序、同一矩阵列序）。"""
    per_probe = {}
    design = list(csv.DictReader(open(outpath('de_queue_design.tsv'), encoding='utf-8'), delimiter='\t'))
    for gse in CHIP_QUEUES:
        probes, mat = m04.load_queue(gse)
        groups_all = m04.load_groups(gse)
        gsm_order = [r['gsm'] for r in design if r['gse'] == gse]
        gmap = {s: groups_all[s] for s in gsm_order}
        beta, t, p, n1, n0 = m04.moderated_t(mat, gmap)
        for i, pb in enumerate(probes):
            per_probe.setdefault(pb, {})[gse] = (float(beta[i]), float(p[i]), n1, n0)
        print(f'  {gse}: {len(probes)} probes (case {n1}/ctl {n0})')
    return per_probe


def combine(per_probe, queues, min_cov=3):
    """04 号合成逻辑的子集版：se 反推 → Stouffer(sqrt n) → DL → 方向一致性。"""
    rows = []
    for pb, d in per_probe.items():
        dd = {g: d[g] for g in queues if g in d}
        if len(dd) < min_cov:
            continue
        gses = sorted(dd)
        betas = [dd[g][0] for g in gses]
        ps = [dd[g][1] for g in gses]
        ns = [dd[g][2] + dd[g][3] for g in gses]
        ses = [abs(b) / max(stats.norm.isf(min(pv, 1 - 1e-12)), 1e-12) if pv > 0 else np.nan
               for b, pv in zip(betas, ps)]
        keep = [i for i, s in enumerate(ses) if not np.isnan(s) and s > 0]
        if len(keep) < min_cov:
            continue
        zs = [betas[i] / ses[i] for i in keep]
        wgt = [np.sqrt(ns[i]) for i in keep]
        z_st = sum(z * w for z, w in zip(zs, wgt)) / np.sqrt(sum(w ** 2 for w in wgt))
        b_dl, se_dl, p_dl, Q, k = m04.dl_random_effect([betas[i] for i in keep], [ses[i] for i in keep])
        pos = sum(1 for b in betas if b > 0)
        consist = max(pos, len(betas) - pos) / len(betas)
        rows.append({'probe': pb, 'n_queues': len(keep), 'beta': b_dl, 'p_DL': p_dl,
                     'consist': consist, 'betas': {gses[i]: betas[i] for i in range(len(gses))}})
    qs = bh([min(r['p_DL'], 1.0) for r in rows])
    for r, q in zip(rows, qs):
        r['q_BH'] = q
        r['tier'] = 'A' if (q < 0.05 and r['consist'] == 1.0) else ('B' if (q < 0.05 and r['consist'] >= 0.8) else ('C' if q < 0.05 else '-'))
    return rows


def sig_same_sign(sub_map, ref_rows):
    """主 meta Tier A 探针在子 meta 中显著且同号的比例。"""
    tot = hit = same = 0
    for r in ref_rows:
        if r['tier'] != 'A':
            continue
        tot += 1
        s = sub_map.get(r['probe'])
        if s and s['q_BH'] < 0.05:
            hit += 1
            if np.sign(s['beta']) == np.sign(r['beta']):
                same += 1
    return {'ref_n': tot, 'sig_in_stratum': hit, 'sig_same_sign': same,
            'retention_pct': round(100 * same / tot, 1) if tot else None}


def main():
    print('per-queue Welch (04-identical):')
    per_probe = compute_per_queue()

    main_rows = combine(per_probe, CHIP_QUEUES)
    main_map = {r['probe']: r for r in main_rows}
    ec_rows = combine(per_probe, ECOLI)
    sa_rows = combine(per_probe, SAUR)
    no_rows = combine(per_probe, NOLIVER)
    ec_map = {r['probe']: r for r in ec_rows}
    sa_map = {r['probe']: r for r in sa_rows}
    no_map = {r['probe']: r for r in no_rows}

    def summ(rows):
        return {'probes_tested': len(rows), 'q05_BH': sum(1 for r in rows if r['q_BH'] < 0.05),
                'tierA': sum(1 for r in rows if r['tier'] == 'A'),
                'tierB': sum(1 for r in rows if r['tier'] == 'B')}

    cross = {
        'main_tierA_in_ecoli': sig_same_sign(ec_map, main_rows),
        'main_tierA_in_saur': sig_same_sign(sa_map, main_rows),
        'main_tierA_in_noliver': sig_same_sign(no_map, main_rows),
    }
    # 病原层之间：两层都显著且同号的共享探针（"共同感染反应轴"的直接检验）
    shared = sorted(set(ec_map) & set(sa_map))
    both_sig = [p for p in shared if ec_map[p]['q_BH'] < 0.05 and sa_map[p]['q_BH'] < 0.05]
    same_sign = [p for p in both_sig if np.sign(ec_map[p]['beta']) == np.sign(sa_map[p]['beta'])]
    # 去肝臂 vs 主 meta 全秩相关
    shared_no = sorted(set(main_map) & set(no_map))
    x = [-np.log10(min(main_map[p]['p_DL'], 1)) for p in shared_no]
    y = [-np.log10(min(no_map[p]['p_DL'], 1)) for p in shared_no]
    rho = stats.spearmanr(x, y).statistic

    # 已发布主 meta 的复算校验（探针层 Tier A 应复现 2,999）
    pub = list(csv.DictReader(open(outpath('de_meta_chip_final.tsv'), encoding='utf-8'), delimiter='\t'))
    pub_tierA = sum(1 for r in pub if r['tier'] == 'A')

    out = {
        'method': 'per-stratum Stouffer(sqrt n) + DerSimonian-Laird, within-stratum coverage >=3 queues, within-stratum BH; identical estimator module as 04_de_meta.py',
        'pathogen_grouping': 'S3 (2026-09-22): E.coli=25319/24217/15025(+75379 Illumina, disjoint probe IDs); S.aureus-related=9685/24560/93082; field=34031/131607 (not probe-combinable: Illumina + RNA-seq)',
        'main_ref_recomputed': summ(main_rows),
        'published_main_tierA': pub_tierA,
        'ecoli': summ(ec_rows),
        'saur': summ(sa_rows),
        'field': {'probes_testable': False,
                  'reason': 'GSE34031 (Illumina) + GSE131607 (RNA-seq): probe/gene ID spaces disjoint, <3 queues; descriptive only'},
        'noliver': summ(no_rows),
        'cross': cross,
        'pathogen_shared': {'shared_probes': len(shared), 'sig_in_both': len(both_sig),
                            'sig_same_sign_both': len(same_sign),
                            'same_sign_pct': round(100 * len(same_sign) / len(both_sig), 1) if both_sig else None},
        'noliver_vs_main_spearman_rho_logp': round(float(rho), 4),
        'note_Danish_pair': 'GSE25319 (liver) and GSE24217 (udder, lobuloalveolar region) are two tissue arms of the same Danish E.coli challenge experiment (same cow IDs in sample metadata); noliver removes the whole arm',
    }
    dump_json(out, 'stratified_meta.json')
    print(json.dumps(out, indent=1))

    # ---- C. Top-30 权重扰动稳定性 ----
    sc = list(csv.DictReader(open(outpath('priority_scoring.tsv'), encoding='utf-8'), delimiter='\t'))
    s1 = np.array([float(r['S1_DE']) for r in sc])
    s3 = np.array([float(r['S3_PIP']) for r in sc])
    s4 = np.array([float(r['S4_network']) for r in sc])
    base = 0.35 * s1 + 0.30 * s3 + 0.10 * s4
    top30 = set(np.argsort(-base)[:30])
    top10 = set(np.argsort(-base)[:10])
    rng = np.random.default_rng(42)
    ov30, ov10, rhos = [], [], []
    for _ in range(1000):
        w = np.array([0.35, 0.30, 0.10]) * rng.uniform(0.5, 1.5, 3)
        w = 0.75 * w / w.sum()
        sc_j = w[0] * s1 + w[1] * s3 + w[2] * s4
        order = np.argsort(-sc_j)
        ov30.append(len(set(order[:30]) & top30))
        ov10.append(len(set(order[:10]) & top10))
        rhos.append(stats.spearmanr(sc_j, base).statistic)

    def stress(zero_idx):
        w = np.array([0.35, 0.30, 0.10], dtype=float)
        w[zero_idx] = 0.0
        w = 0.75 * w / w.sum()
        order = np.argsort(-(w[0] * s1 + w[1] * s3 + w[2] * s4))
        return {'top30_overlap': len(set(order[:30]) & top30), 'top10_overlap': len(set(order[:10]) & top10)}

    wp = {
        'draws': 1000, 'seed': 42,
        'perturbation': 'S1/S3/S4 weights independently x U(0.5,1.5), renormalized to 0.75 (S2 constant, not rank-active)',
        'top30_overlap': {'median': float(np.median(ov30)), 'min': int(np.min(ov30)), 'max': int(np.max(ov30))},
        'top10_overlap': {'median': float(np.median(ov10)), 'min': int(np.min(ov10)), 'max': int(np.max(ov10))},
        'rank_spearman_median': round(float(np.median(rhos)), 4),
        'stress': {'no_DE(w1=0)': stress(0), 'no_PIP(w3=0)': stress(1), 'no_network(w4=0)': stress(2)},
        'note': 'ranking base excludes the constant S2 term (0.25*0.5 for every scored gene)',
    }
    dump_json(wp, 'weight_perturbation.json')
    print(json.dumps(wp, indent=1))


if __name__ == '__main__':
    main()
