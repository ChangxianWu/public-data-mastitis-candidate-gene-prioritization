#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 1d: BH 校正 + 方向一致性分层汇总（DE meta 芯片轨收口表）。

输出 work/output/de_meta_chip_final.tsv：
  probe, n_queues, beta_meta, p_DL, q_BH, sign_consistency, tier
  tier A = q<0.05 且一致性 1.0（5/5 或全队列同号）
  tier B = q<0.05 且一致性 >= 0.8
  tier C = 其余 q<0.05
"""
import csv, json, os, sys
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/output'


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    q = [1.0] * n
    running = 1.0
    for rank, i in enumerate(reversed(order), 1):
        running = min(running, ps[i] * n / (n - rank + 1))
        q[i] = running
    return q


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, 'de_meta_chip_summary.tsv'), encoding='utf-8'), delimiter='\t'))
    ps = [min(float(r['p_DL']), 1.0) for r in rows]
    qs = bh(ps)
    tier_a = tier_b = tier_c = 0
    with open(os.path.join(OUT, 'de_meta_chip_final.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['probe', 'n_queues', 'beta_meta', 'p_DL', 'q_BH', 'sign_consistency', 'tier'])
        for r, q in zip(rows, qs):
            consist = float(r['sign_consistency'])
            if q < 0.05 and consist == 1.0:
                tier, tier_a = 'A', tier_a + 1
            elif q < 0.05 and consist >= 0.8:
                tier, tier_b = 'B', tier_b + 1
            elif q < 0.05:
                tier, tier_c = 'C', tier_c + 1
            else:
                tier = '-'
            w.writerow([r['probe'], r['n_queues'], r['beta_meta'], r['p_DL'],
                        f'{q:.4g}', r['sign_consistency'], tier])
    summary = {'probes_meta': len(rows), 'q_BH_0.05': tier_a + tier_b + tier_c,
               'tierA_q05_consist100': tier_a, 'tierB_q05_consist80': tier_b, 'tierC': tier_c}
    json.dump(summary, open(os.path.join(OUT, 'de_meta_chip_final_summary.json'), 'w', encoding='utf-8'), indent=1)
    print(summary)


if __name__ == '__main__':
    main()
