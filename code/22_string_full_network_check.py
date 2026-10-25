#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 4a-check: STRING v11.5 官方物种全网络核对（评审 M8c）。

主网络层仍用冻结的 API 快照（string_edges_raw.json，随包分发）以避免端点不稳定；
本脚本从 STRING 官方下载站取 9913 全网络 + 蛋白注释，核验：
  1. 文件 SHA-256（随包记录）
  2. 11 个网络正分基因在官方全网络中的物种全网络加权度（score >= 400）
  3. 与 API 子网度的秩一致性（报告两者，明确二者口径不同）

输出 work/output/string_full_network_check.json
运行：python work/code/22_string_full_network_check.py
"""
import csv, gzip, hashlib, json, os, sys
from collections import defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.normpath(os.path.join(HERE, 'output'))
DATA = os.path.join(HERE, 'data', 'string')
LINKS = os.path.join(DATA, '9913.protein.links.v11.5.txt.gz')
INFO = os.path.join(DATA, '9913.protein.info.v11.5.txt.gz')


def outpath(name):
    base = os.path.normpath(OUT)
    p = os.path.normpath(base + os.sep + name)
    if not p.startswith(base + os.sep):
        raise SystemExit('output path escaped OUT dir: ' + name)
    return p


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 20), b''):
            h.update(ch)
    return h.hexdigest()


def main():
    if not (os.path.isfile(LINKS) and os.path.isfile(INFO)):
        raise SystemExit('missing STRING download files under work/data/string')
    # protein id -> preferred name
    p2n = {}
    with gzip.open(INFO, 'rt', encoding='utf-8', errors='replace') as fh:
        hdr = fh.readline()
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2 and p[1]:
                p2n[p[0]] = p[1]
    name2prot = defaultdict(list)
    for pid, nm in p2n.items():
        name2prot[nm].append(pid)
    # seeds = genes with S4>0 in the published scoring
    seeds = []
    for r in csv.DictReader(open(outpath('priority_scoring.tsv'), encoding='utf-8'), delimiter='\t'):
        if float(r['S4_network']) > 0 and r['symbol']:
            seeds.append((r['symbol'], int(r['STRING_degree'])))
    seeds = sorted(set(seeds))
    want = {s for s, _ in seeds}
    wdeg = defaultdict(float)
    deg = defaultdict(int)
    with gzip.open(LINKS, 'rt', encoding='utf-8', errors='replace') as fh:
        hdr = fh.readline()
        for line in fh:
            p = line.rstrip('\n').split()
            if len(p) < 3:
                continue
            a, b, s = p[0], p[1], float(p[2])
            if s < 400:
                continue
            na, nb = p2n.get(a), p2n.get(b)
            if na in want:
                wdeg[na] += s
                deg[na] += 1
            if nb in want:
                wdeg[nb] += s
                deg[nb] += 1
    rows = []
    for sym, api_deg in seeds:
        rows.append({'symbol': sym, 'api_subnetwork_degree': api_deg,
                     'full_network_degree': deg.get(sym, 0),
                     'full_network_weighted_degree': round(wdeg.get(sym, 0.0), 1),
                     'protein_ids': name2prot.get(sym, [])[:3]})
    rank_api = {r['symbol']: i + 1 for i, r in enumerate(sorted(rows, key=lambda x: -x['api_subnetwork_degree']))}
    rank_full = {r['symbol']: i + 1 for i, r in enumerate(sorted(rows, key=lambda x: -x['full_network_weighted_degree']))}
    same_order = all(rank_api[s] == rank_full[s] for s in rank_api)
    out = {
        'files': {
            os.path.basename(LINKS): {'sha256': sha256(LINKS), 'bytes': os.path.getsize(LINKS)},
            os.path.basename(INFO): {'sha256': sha256(INFO), 'bytes': os.path.getsize(INFO)},
        },
        'score_threshold': 400,
        'seeds': rows,
        'rank_same_order_api_vs_full': same_order,
        'note': 'API subnetwork degree counts edges among the returned seed subnetwork; full-network degree counts all species-wide partners (score>=400). Different quantities; only rank order is compared.',
    }
    Path(outpath('string_full_network_check.json')).write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(out, indent=1, ensure_ascii=False)[:2500])


if __name__ == '__main__':
    main()
