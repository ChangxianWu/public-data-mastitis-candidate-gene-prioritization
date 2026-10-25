#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 3a: 基因锚点表生成（复核修复轮补齐的生成脚本；此前 de_meta_gene_anchored.tsv 无落盘脚本）。

规则（由 2026-09-19 产物逆向，逐字段对拍 0 差异确认）：
  1. GPL2112_probe2gene_final.tsv（platform_direct + RefSeq efetch 两通道合并表）
     提供探针 → (gene_symbol, entrez_id) 映射。
  2. **按 gene_symbol 归组**（非 Entrez：同一 Ensembl 基因经不同通道可获不同 Entrez，
     按符号归组是本表历史口径）；空符号探针（5 个）不进入基因层。
  3. 每符号取 p_DL 最小的探针为锚（β/p/一致性/n_queues 取锚探针），
     n_probes_gene = 该符号映射探针数，n_same_direction = 与锚探针同号的探针数。
  4. 按 p_DL 升序输出（1,886 基因）。

输入：work/data/platform/GPL2112_probe2gene_final.tsv、work/output/de_meta_chip_final.tsv
输出：work/output/de_meta_gene_anchored.tsv
运行方式：在 VM2 课题根目录执行  python work/code/07b_gene_anchor.py
"""
import csv, os, sys
from collections import defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
ROOT = Path('.').resolve()
PLAT = ROOT / 'work' / 'data' / 'platform' / 'GPL2112_probe2gene_final.tsv'
OUTDIR = ROOT / 'work' / 'output'
for p in (PLAT, OUTDIR):
    if not p.exists():
        sys.exit(f'请在 VM2 课题根目录运行（缺少 {p}）')


def main():
    p2g = {}
    for r in csv.DictReader(open(PLAT, encoding='utf-8'), delimiter='\t'):
        p2g[r['probe']] = (r['gene_symbol'], r['entrez_id'])

    chip = {r['probe']: r for r in csv.DictReader(open(OUTDIR / 'de_meta_chip_final.tsv',
                                                      encoding='utf-8'), delimiter='\t')}
    by_sym = defaultdict(list)
    n_blank = 0
    for pb, (sym, ent) in p2g.items():
        if pb not in chip:
            continue
        if not sym.strip():
            n_blank += 1
            continue
        by_sym[sym].append((pb, ent))
    print(f'symbols: {len(by_sym)}  probes skipped (blank symbol): {n_blank}')

    rows = []
    for sym, lst in by_sym.items():
        lst.sort(key=lambda x: float(chip[x[0]]['p_DL']))
        pb, ent = lst[0]
        b_anchor = float(chip[pb]['beta_meta'])
        rows.append({'gene_symbol': sym, 'entrez_id': ent, 'probe': pb,
                     'n_queues': chip[pb]['n_queues'], 'beta_meta': chip[pb]['beta_meta'],
                     'p_DL': chip[pb]['p_DL'], 'sign_consistency': chip[pb]['sign_consistency'],
                     'n_probes_gene': len(lst),
                     'n_same_direction': sum(1 for p, _ in lst
                                             if float(chip[p]['beta_meta']) * b_anchor > 0)})
    # 稳定排序：p_DL 升序；同 p 的并列按平台映射文件的符号出现顺序（与历史产物一致）
    rows.sort(key=lambda r: float(r['p_DL']))

    out = OUTDIR / 'de_meta_gene_anchored.tsv'
    with open(out, 'w', encoding='utf-8', newline='') as f:
        wtr = csv.DictWriter(f, fieldnames=['gene_symbol', 'entrez_id', 'probe', 'n_queues',
                                            'beta_meta', 'p_DL', 'sign_consistency',
                                            'n_probes_gene', 'n_same_direction'], delimiter='\t')
        wtr.writeheader()
        wtr.writerows(rows)
    up = sum(1 for r in rows if float(r['beta_meta']) > 0)
    print(f'wrote {out.name}: {len(rows)} genes (up {up} / down {len(rows) - up})')


if __name__ == '__main__':
    main()
