#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 GEO 下载器（DE meta 层 Stage 1）。

关键工程决定（2026-09-17 解冻依据）：NCBI 域经本机代理（127.0.0.1:10808）TLS 握手必超时，
直连 4.5s 即 200——本下载器对 NCBI 域**强制绕过代理**（ProxyHandler({})），对其他域走环境默认。

功能：
  1. 代表集 GSE 的 series_matrix 下载（小文件，提取元数据 + 表达矩阵定位）
  2. 元数据汇总表（GSE × 平台 × 样本数 × 组织 × 病原分组字段）
用法: python work/code/01_geo_download.py
"""
import json, os, re, sys, time, urllib.request, gzip, io

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # work/
OUT = os.path.join(HERE, 'output')
DATA = os.path.join(HERE, 'data', 'geo')
os.makedirs(DATA, exist_ok=True)

NCBI_OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
NCBI_OPENER.addheaders = [('User-Agent', 'VM2-gate/1.0 (DE-meta stage1)')]

REPRESENTATIVES = ['GSE202849', 'GSE256521', 'GSE334132', 'GSE334081', 'GSE296484']


def gse_url(gse, fn):
    # GEO 惯例：目录 = GSE + 前三位数字 + nnn（GSE202849 -> GSE202nnn）
    n = int(re.match(r'GSE(\d+)', gse).group(1))
    nn = 'GSE{}nnn'.format(n // 1000)
    return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{nn}/{gse}/matrix/{fn}'


def fetch(url, tries=4, timeout=90):
    last = None
    for i in range(tries):
        try:
            return NCBI_OPENER.open(url, timeout=timeout).read()
        except Exception as e:
            last = e
            time.sleep(3 + i * 4)
    raise last


def main():
    meta = {}
    for gse in REPRESENTATIVES:
        fn = f'{gse}_series_matrix.txt.gz'
        dst = os.path.join(DATA, fn)
        if not os.path.exists(dst) or os.path.getsize(dst) < 500:
            url = gse_url(gse, fn)
            print(gse, 'downloading...')
            b = fetch(url)
            open(dst, 'wb').write(b)
            time.sleep(1.2)
        raw = gzip.open(dst, 'rt', encoding='utf-8', errors='replace').read()
        rec = {'gse': gse, 'bytes': os.path.getsize(dst)}
        for key in ('!Series_title', '!Series_type', '!Series_platform_id', '!Series_sample_taxid',
                    '!Series_overall_design'):
            m = re.search(rf'^{re.escape(key)}\s*=\s*(.*)$', raw, re.M)
            if m:
                rec[key] = m.group(1).strip().strip('"')[:200]
        # 样本数与样本表头
        ids = re.search(r'!Sample_geo_accession\s*=\s*(.*)', raw)
        rec['n_samples'] = len(ids.group(1).split('\t')) if ids else 0
        chars = re.findall(r'!Sample_characteristics_ch1\s*=\s*(.*)', raw)
        rec['n_char_fields'] = len(chars)
        rec['char_first_tokens'] = sorted({c.split(':')[0].strip() for c in chars if ':' in c})[:8]
        # 表达矩阵存在性
        rec['has_series_matrix_expr'] = 'series_matrix_table_begin' in raw
        # 平台类型判定
        plat = rec.get('!Series_platform_id', '')
        rec['is_rnaseq'] = bool(re.search(r'GPL\d+', plat)) and ('High throughput sequencing' in rec.get('!Series_type', '') or 'RNA' in rec.get('!Series_type', ''))
        meta[gse] = rec
        print(f"  {gse}: samples={rec['n_samples']} chars={rec['char_first_tokens']} expr_in_matrix={rec['has_series_matrix_expr']}")
        time.sleep(1.0)
    out = os.path.join(OUT, 'representative_series_matrix_summary.json')
    json.dump(meta, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('saved', out)


if __name__ == '__main__':
    main()
