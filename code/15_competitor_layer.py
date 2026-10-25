#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 6b: 竞文对比层（文献基因集）提取（复核修复轮补齐的生成脚本）。

竞文 = Cai Z et al. 2024（PMID 39009986）六品种基因级 meta 的补充表
MOESM1（12711_2024_920_MOESM1_ESM.docx）中的四张基因表 T5/T6/T7/T8（docx 内序）：
  T5  46 行  gene 列 = 符号（SOX14 ...）
  T6  62 行  GENE ID 列 = Ensembl ENSBTAG（RBM48 ...）
  T7   7 行  同 T6（MTAG_CM）
  T8  65 行  GENE ID 列 = Ensembl ENSBTAG
提取规则：每表取基因标识列非空且不为 'NA' 的值；T5 用符号原样、T6/T7/T8 用
ENSBTAG ID 原样；去重合并 → 文献对比层（148 基因）；与本项目高置信层
（de_eqtl_pip_integrated.tsv 中 PIP>=0.5 & p<0.05 的 Ensembl ID）求交。

输入：work/data_raw/competitor_39009986/12711_2024_920_MOESM1_ESM.docx
      work/output/de_eqtl_pip_integrated.tsv
输出：work/output/mini_qtl_gene_baseline.json / mini_qtl_intersection.json
运行方式：在 VM2 课题根目录执行  python work/code/15_competitor_layer.py
"""
import json, sys
from pathlib import Path

import docx

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
ROOT = Path('.').resolve()
SRC = ROOT / 'work' / 'data_raw' / 'competitor_39009986' / '12711_2024_920_MOESM1_ESM.docx'
OUTDIR = ROOT / 'work' / 'output'
if not SRC.is_file() or not OUTDIR.is_dir():
    sys.exit('请在 VM2 课题根目录运行（需 work/data_raw/competitor_39009986 与 work/output）')


def cell(r, i):
    try:
        v = r.cells[i].text.strip()
    except Exception:
        v = ''
    return '' if v.upper() == 'NA' else v


def gene_token(v):
    """基因标识归一：ENSBTAG 原样；符号统一大写；剔除表头残留。"""
    v = v.strip()
    if not v or v.upper() in ('NA', 'GENE', 'GENE ID', 'GENE SYMBOL'):
        return None
    return v if v.startswith('ENSBTAG') else v.upper()


def main():
    d = docx.Document(str(SRC))
    baseline = {}
    # T5 只有 gene 符号列；T6/T7/T8 同时取 GENE ID（Ensembl）与 GENE（符号）两列
    plan = ((5, (1,)), (6, (1, 2)), (7, (1, 2)), (8, (1, 2)))
    for ti, cols in plan:
        t = d.tables[ti]
        n = 0
        for r in t.rows[1:]:
            for col in cols:
                g = gene_token(cell(r, col))
                if not g:
                    continue
                baseline.setdefault(g, {'sources': []})
                if f'T{ti}' not in baseline[g]['sources']:
                    baseline[g]['sources'].append(f'T{ti}')
                n += 1
        print(f'T{ti}: {n} 个基因标识 token')

    uniq = sorted(baseline)
    print('unique competitor genes:', len(uniq))

    lines = (OUTDIR / 'de_eqtl_pip_integrated.tsv').read_text(encoding='utf-8').splitlines()[1:]
    hi = set(l.split('\t')[0] for l in lines
             if float(l.split('\t')[6]) >= 0.5 and float(l.split('\t')[4]) < 0.05)
    inter = sorted(hi & set(uniq))
    print('our high-confidence layer:', len(hi), ' intersection:', inter)

    (OUTDIR / 'mini_qtl_gene_baseline.json').write_text(
        json.dumps(baseline, indent=1, ensure_ascii=False), encoding='utf-8')
    (OUTDIR / 'mini_qtl_intersection.json').write_text(
        json.dumps(inter, indent=1), encoding='utf-8')
    print('wrote mini_qtl_gene_baseline.json / mini_qtl_intersection.json')


if __name__ == '__main__':
    main()
