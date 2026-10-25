#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 2a: 平台注释下载 + 探针→基因映射（GPL2112/6082/11648）。

GPL*_family.soft.gz 的 platform_table 含 ID / GeneSymbol / GeneID / GenBank 等列。
输出 work/data/platform/<GPL>_probe2gene.tsv (probe, gene_symbol, entrez_id)
"""
import gzip, io, os, re, sys, time, urllib.request

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDIR = os.path.join(HERE, 'data', 'platform')
os.makedirs(PDIR, exist_ok=True)

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
OPENER.addheaders = [('User-Agent', 'VM2-gate/1.0')]

GPL = {'GPL2112': 'GPL2nnn', 'GPL6082': 'GPL6nnn', 'GPL11648': 'GPL11nnn'}


def fetch(url, tries=4, timeout=120):
    last = None
    for i in range(tries):
        try:
            r = OPENER.open(url, timeout=timeout)
            chunks = []
            while True:
                c = r.read(1 << 18)
                if not c:
                    break
                chunks.append(c)
            return b''.join(chunks)
        except Exception as e:
            last = e
            print('  retry', i + 1, type(e).__name__, flush=True)
            time.sleep(4 + i * 4)
    raise last


def parse_soft(text):
    """从 family.soft 提取 platform_table：probe -> (symbol, entrez)"""
    mapping = {}
    cols = None
    in_table = False
    for ln in text.split('\n'):
        if ln.startswith('^PLATFORM'):
            in_table = False
            cols = None
        elif ln.startswith('!platform_table_begin'):
            in_table = True
        elif ln.startswith('!platform_table_end'):
            in_table = False
        elif in_table:
            if cols is None:
                cols = [c.strip() for c in ln.lstrip('#!').split('\t')]
                continue
            parts = ln.rstrip('\r').split('\t')
            if len(parts) < 2:
                continue
            row = dict(zip(cols, parts))
            pid = row.get('ID', '').strip('"')
            sym = (row.get('GeneSymbol') or row.get('Gene Symbol') or row.get('GENE_SYMBOL') or '').strip().strip('"')
            entrez = (row.get('GeneID') or row.get('Entrez_Gene_ID') or row.get('ENTREZ_GENE_ID') or '').strip().strip('"')
            if pid:
                mapping[pid] = (sym, entrez)
    return mapping


def main():
    for gpl, nn in GPL.items():
        dst = os.path.join(PDIR, f'{gpl}_probe2gene.tsv')
        if os.path.exists(dst) and os.path.getsize(dst) > 1000:
            print(gpl, 'already mapped')
            continue
        url = f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{nn}/{gpl}/soft/{gpl}_family.soft.gz'
        print(gpl, 'downloading family.soft...', flush=True)
        b = fetch(url)
        text = gzip.decompress(b).decode('utf-8', errors='replace')
        mp = parse_soft(text)
        n_sym = sum(1 for s, _ in mp.values() if s)
        with open(dst, 'w', encoding='utf-8', newline='') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['probe', 'gene_symbol', 'entrez_id'])
            for pid, (sym, entrez) in sorted(mp.items()):
                w.writerow([pid, sym, entrez])
        print(f'  {gpl}: {len(mp)} probes, {n_sym} with gene symbol')


import csv
if __name__ == '__main__':
    main()
