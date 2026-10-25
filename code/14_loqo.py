#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 4b: LOQO（leave-one-queue-out）稳定性复算（复核修复轮补齐的生成脚本）。

对 Tier A（q<0.05 且全队列同号，2,999 探针）逐个剔除一支芯片队列后按 04 号同款
流程重算 DL 随机效应，统计"仍显著且方向不变"的存活数。

2026-10-25 更正：仅对真正进入探针级 meta 的五个 GPL2112 队列做剔除。GSE9685 与
GSE34031 的沉积矩阵无探针标识（数字行号），从未进入 meta，对其"剔除"是空操作
（旧版因此产出 99.1% 两档虚高），已从队列表中移除；剔除后维持 >=3 队列覆盖约束。

输入：work/output/de_meta_chip_per_queue.tsv（每队列 β/t/p）、de_meta_chip_final.tsv（Tier A 名单）
输出：work/output/loqo_stability.json（5 队列 × {tested, survived, rate}）
运行方式：在 VM2 课题根目录执行  python work/code/14_loqo.py
"""
import csv, json, os, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
OUTDIR = Path('work/output').resolve()
if not OUTDIR.is_dir():
    sys.exit('请在 VM2 课题根目录运行')


def rd(name):
    p = OUTDIR / name
    if not p.is_file():
        sys.exit(f'缺少输入: {p}')
    return p


def dl_random_effect(betas, ses):
    betas = np.asarray(betas, float)
    ses = np.asarray(ses, float)
    ok = ~np.isnan(betas) & ~np.isnan(ses)
    betas, ses = betas[ok], ses[ok]
    k = len(betas)
    if k == 0:
        return np.nan, np.nan, 0
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
    return b_re, p, k


def main():
    per_q = defaultdict(dict)   # probe -> {gse: (beta, p)}
    for r in csv.DictReader(open(rd('de_meta_chip_per_queue.tsv'),
                                 encoding='utf-8'), delimiter='\t'):
        if r['p'] not in ('NA', ''):
            per_q[r['probe']][r['gse']] = (float(r['beta']), float(r['p']))
    tier_a = [(r['probe'], float(r['beta_meta']))
              for r in csv.DictReader(open(rd('de_meta_chip_final.tsv'),
                                           encoding='utf-8'), delimiter='\t') if r['tier'] == 'A']
    print('Tier A probes:', len(tier_a))

    queues = ['GSE24560', 'GSE25319', 'GSE24217', 'GSE93082', 'GSE15025']
    res = {}
    for drop in queues:
        survived_sig_dir = survived_sig = 0
        tested = 0
        for pb, b_full in tier_a:
            d = per_q.get(pb, {})
            items = [(g, b, p) for g, (b, p) in d.items() if g != drop]
            if len(items) < 3:
                continue
            betas = [b for _, b, _ in items]
            ses = []
            for _, b, p in items:
                z = stats.norm.isf(min(max(p, 1e-12), 1 - 1e-12))
                ses.append(abs(b) / z if z > 0 else np.nan)
            keep = [i for i, s in enumerate(ses) if not np.isnan(s) and s > 0]
            if len(keep) < 2:
                continue
            tested += 1
            b_re, p_re, _ = dl_random_effect([betas[i] for i in keep], [ses[i] for i in keep])
            if p_re is not None and not np.isnan(p_re) and p_re < 0.05:
                survived_sig += 1
                if not np.isnan(b_re) and b_re * b_full > 0:
                    survived_sig_dir += 1
        res[drop] = {'tested': tested, 'survived': survived_sig_dir,
                     'rate': round(survived_sig_dir / tested, 3) if tested else None}
        print(f'{drop}: tested={tested} survived(sig&dir)={survived_sig_dir} '
              f'(sig-only={survived_sig}) rate={res[drop]["rate"]}')

    (OUTDIR / 'loqo_stability.json').write_text(json.dumps(res, indent=1), encoding='utf-8')
    print('wrote loqo_stability.json')


if __name__ == '__main__':
    main()
