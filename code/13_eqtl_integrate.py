#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 3b: 乳腺 cis-eQTL / fine-mapping 整合层（复核修复轮补齐的生成脚本）。

背景（2026-10-06 复核）：本层产物原先无落盘脚本。本脚本以两种口径重建，并作为
de_eqtl_overlap_* / de_eqtl_pip_* / de_eqtl_pip_integrated.tsv 的唯一生成者：

  strict （默认，稿件口径）  integrated = DE合并基因 ∩ eQTL显著集 ∩ finemapping
                             = 3,395 基因（与摘要"3,395 with fine-mapped lead variants"
                               与高置信层定义 DE×eQTL presence×PIP 完全一致）
  extended（legacy 对拍用）  integrated = DE合并基因 ∩ finemapping
                             = 3,404 基因（2026-09-19 手工版口径；9 个基因无显著
                               eQTL 证据，其中 2 个曾进入 288 层——复核报告 M2，
                               2026-10-07 起弃用，此口径仅留作复现对拍）

输入：
  work/output/de_meta_gene_combined.tsv              （08 号产物，24,597 基因）
  $CATTLE_GTEX_MAMMARY_DIR/Mammary_Gland_significant_eQTL.txt.gz
  $CATTLE_GTEX_MAMMARY_DIR/Mammary_Gland_finemapping_eQTL.txt.gz
  环境变量 CATTLE_GTEX_MAMMARY_DIR 指向 CattleGTEx Mammary_Gland 快照目录
  （本工作区的快照在 VM1/work/data_raw/cattlegtex_mammary_eqtl/ 下；获取方法见
   VM1 S12 数据访问手册。注意 significant 文件名为 .gz 实为纯 TSV——分发陷阱）。
输出：
  work/output/de_eqtl_overlap_genes.txt / _summary.json / de_eqtl_pip_summary.json
  work/output/de_eqtl_pip_integrated.tsv（--caliber 指定口径）
"""
import csv, gzip, json, os, sys
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, 'output')

SRC = os.environ.get('CATTLE_GTEX_MAMMARY_DIR', '')
if not SRC:
    sys.exit('需要环境变量 CATTLE_GTEX_MAMMARY_DIR（CattleGTEx Mammary_Gland 快照目录，见 VM1 S12 手册）')
SRC = os.path.realpath(SRC)
if not os.path.isdir(SRC):
    sys.exit(f'快照目录不存在: {SRC}')
SIG = os.path.join(SRC, 'Mammary_Gland_significant_eQTL.txt.gz')
FM = os.path.join(SRC, 'Mammary_Gland_finemapping_eQTL.txt.gz')
for p in (SIG, FM):
    if not os.path.isfile(p):
        sys.exit(f'缺少输入文件: {p}')


def opener(p):
    with open(p, 'rb') as f:
        magic = f.read(2)
    return gzip.open if magic == b'\x1f\x8b' else open


def rows_of(path):
    op = opener(path)
    kw = {} if op is gzip.open else {'encoding': 'utf-8', 'errors': 'replace'}
    with op(path, 'rt', **kw) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            yield line.rstrip('\n').split('\t')


def main():
    caliber = 'strict'
    argv = sys.argv[1:]
    if '--caliber' in argv:
        caliber = argv[argv.index('--caliber') + 1]
    assert caliber in ('strict', 'extended'), caliber

    # 1) DE 合并基因（08 号产物）
    de = list(csv.DictReader(open(os.path.join(OUT, 'de_meta_gene_combined.tsv'),
                                  encoding='utf-8'), delimiter='\t'))
    print('DE combined genes:', len(de))

    # 2) eQTL 显著集（7,206）
    sig = set()
    for c in rows_of(SIG):
        if c:
            sig.add(c[0])
    print('significant eQTL genes:', len(sig))

    # 3) DE ∩ 显著集 → overlap 明细（列序与历史版一致：gene,symbol,beta,p,nq,sources）
    overlap = [r for r in de if r['gene'] in sig]
    overlap.sort(key=lambda r: float(r['p_DL']) if r['p_DL'] else 2)
    with open(os.path.join(OUT, 'de_eqtl_overlap_genes.txt'), 'w', encoding='utf-8', newline='') as f:
        for r in overlap:
            f.write('\t'.join([r['gene'], r['symbol'], r['beta_meta'], r['p_DL'],
                               r['n_queues'], r['sources']]) + '\n')
    json.dump({'eqtl_genes': len(sig), 'de_genes': len(de), 'overlap': len(overlap)},
              open(os.path.join(OUT, 'de_eqtl_overlap_summary.json'), 'w', encoding='utf-8'), indent=1)
    print('DE ∩ significant:', len(overlap))

    # 4) finemapping：每基因 max PIP + lead SNP
    maxpip, lead = {}, {}
    n_fm = 0
    for c in rows_of(FM):
        if len(c) < 5:
            continue
        n_fm += 1
        g, snp = c[1], c[2]
        try:
            v = float(c[4])
        except ValueError:
            continue
        if g not in maxpip or v > maxpip[g]:
            maxpip[g] = v
            lead[g] = snp
    print(f'finemapping rows: {n_fm}, genes: {len(maxpip)}')

    # 5) integrated 表（口径开关）
    if caliber == 'strict':
        integ = [r for r in overlap if r['gene'] in maxpip]
    else:
        integ = [r for r in de if r['gene'] in maxpip]
    integ.sort(key=lambda r: float(r['p_DL']) if r['p_DL'] else 2)
    out_name = 'de_eqtl_pip_integrated.tsv' if caliber == 'strict' \
        else '_integrated_extended_check.tsv'
    with open(os.path.join(OUT, out_name), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['gene', 'symbol', 'n_queues', 'beta_meta', 'p_DL', 'sources', 'max_PIP', 'lead_SNP'])
        for r in integ:
            w.writerow([r['gene'], r['symbol'], r['n_queues'], f"{float(r['beta_meta']):.4g}",
                        f"{float(r['p_DL']):.4g}", r['sources'],
                        repr(maxpip[r['gene']]), lead[r['gene']]])

    pip_cov = sum(1 for r in overlap if r['gene'] in maxpip)
    top10 = sorted(maxpip.items(), key=lambda kv: -kv[1])[:10]
    json.dump({'overlap_genes_with_pip': pip_cov,
               'max_pip_top10': [[g, v] for g, v in top10]},
              open(os.path.join(OUT, 'de_eqtl_pip_summary.json'), 'w', encoding='utf-8'), indent=1)
    hi = sum(1 for r in integ if maxpip[r['gene']] >= 0.5 and float(r['p_DL']) < 0.05)
    strict9 = sum(1 for r in integ if maxpip[r['gene']] >= 0.9 and float(r['p_DL']) < 0.05)
    print(f'[{caliber}] integrated rows: {len(integ)}  high-conf(PIP>=.5&p<.05): {hi}  strict(>=.9&p<.05): {strict9}')
    print('wrote:', out_name)


if __name__ == '__main__':
    main()
