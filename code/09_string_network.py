#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 4a: STRING 网络层（高置信层种子基因，strict 口径 = PIP>=0.5 & p<0.05）。

- identifiers 用 gene symbol（STRING 牛库识别 symbol；无符号基因回落 Ensembl gene ID）
- 每批 60 基因（API 稳妥批）；network 端点取全互作对 → 度中心性 + 分数
- STRING 证书过期 → unverified SSL 上下文（管线注意点）

可复现性（2026-10-07 复核修复）：STRING network 端点每次返回的互作集合不保证同一
（同一批种子两次调用可得到不同边集：实测 9,984 vs 10,962），因此本层**默认读缓存**
work/output/string_edges_raw.json（首次抓取的冻结快照，随包分发），只有加 --refresh
才重新调 API 并覆盖缓存。稿件报的边/节点数即冻结快照口径。

输出 work/output/string_network.tsv + string_degree.tsv
用法：python work/code/09_string_network.py [--refresh]
"""
import csv, json, os, sys, time, urllib.request, urllib.parse, ssl
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, 'output')
CACHE = os.path.join(OUT, 'string_edges_raw.json')

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}), urllib.request.HTTPSHandler(context=CTX))
OPENER.addheaders = [('User-Agent', 'VM2-gate/1.0')]


def fetch_edges(ids):
    B = 60
    edges = []
    for i in range(0, len(ids), B):
        chunk = ids[i:i + B]
        params = {'identifiers': '%0d'.join(chunk), 'species': 9913, 'limit': 200,
                  'caller_identity': 'VM2_mastitis'}
        url = 'https://string-db.org/api/tsv-no-header/network?' + urllib.parse.urlencode(params)
        for a in range(4):
            try:
                t = OPENER.open(url, timeout=120).read().decode()
                break
            except Exception as e:
                print('  retry', a + 1, type(e).__name__, flush=True)
                time.sleep(5 + a * 5)
        else:
            print(f'  batch {i//B}: FAILED, skipping')
            continue
        for ln in t.strip().split('\n'):
            p = ln.split('\t')
            if len(p) >= 10:
                edges.append((p[2], p[3], float(p[5])))
        print(f'  batch {i//B + 1}/{(len(ids)-1)//B + 1}: cum edges {len(edges)}', flush=True)
        time.sleep(1.5)
    return edges


def main():
    refresh = '--refresh' in sys.argv
    ids = []
    for r in csv.DictReader(open(os.path.join(OUT, 'de_eqtl_pip_integrated.tsv'), encoding='utf-8'), delimiter='\t'):
        if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05:
            ids.append(r['symbol'] or r['gene'])
    ids = sorted(set(ids))
    print('network seed genes:', len(ids))

    if not refresh and os.path.exists(CACHE):
        edges = [tuple(x) for x in json.load(open(CACHE))]
        print(f'using frozen edge cache ({os.path.basename(CACHE)}, {len(edges)} raw edges)')
    else:
        edges = fetch_edges(ids)
        json.dump(edges, open(CACHE, 'w'))
        print(f'fetched and cached {len(edges)} raw edges -> {os.path.basename(CACHE)}')

    uniq = {(min(a, b), max(a, b), s) for a, b, s in edges}
    deg = defaultdict(int)
    wdeg = defaultdict(float)
    for a, b, s in uniq:
        deg[a] += 1
        deg[b] += 1
        wdeg[a] += s
        wdeg[b] += s
    with open(os.path.join(OUT, 'string_network.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['node_a', 'node_b', 'score'])
        for a, b, s in sorted(uniq, key=lambda x: -x[2]):
            w.writerow([a, b, s])
    with open(os.path.join(OUT, 'string_degree.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['gene', 'degree', 'weighted_degree'])
        for g, d in sorted(deg.items(), key=lambda x: -x[1]):
            w.writerow([g, d, round(wdeg[g], 3)])
    print(f'edges: {len(uniq)} unique; nodes with degree: {len(deg)}')


if __name__ == '__main__':
    main()
