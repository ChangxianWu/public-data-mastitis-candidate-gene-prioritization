#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Resolve official gene symbols for the integrated gene set (Ensembl REST).

Review round 2 (A2): the integrated/priority layers carried only platform-derived symbols, so
genes without a platform symbol were presented as "unnamed" even when Ensembl holds a display
name. This stage resolves every integrated Ensembl gene ID against the Ensembl REST
`lookup/id` endpoint (batch POST, species bos_taurus, ARS-UCD2.0) and records symbol, biotype
and description. The mapping is used for display (Tables 3-4, Figure 6) and for reporting the
resolution rate; it does not change any scoring input.

Output: work/output/ensembl_symbol_resolution.json
Run: python work/code/23_symbol_resolution.py
"""
import csv, json, os, ssl, sys, time, urllib.request
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = Path(__file__).resolve().parents[1]
OUT = HERE / 'output'
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
OP = urllib.request.build_opener(urllib.request.ProxyHandler({}), urllib.request.HTTPSHandler(context=CTX))


def lookup(ids):
    body = json.dumps({'ids': ids}).encode()
    req = urllib.request.Request(
        'https://rest.ensembl.org/lookup/id?species=bos_taurus', data=body,
        headers={'Content-Type': 'application/json', 'Accept': 'application/json',
                 'User-Agent': 'VM2/1.0'})
    for a in range(4):
        try:
            return json.loads(OP.open(req, timeout=120).read())
        except Exception as e:
            print('  retry', a + 1, type(e).__name__)
            time.sleep(3 * (a + 1))
    raise SystemExit('Ensembl lookup failed')


def main():
    rows = list(csv.DictReader(open(OUT / 'de_eqtl_pip_integrated.tsv', encoding='utf-8'), delimiter='\t'))
    ids = sorted({r['gene'] for r in rows})
    res = {}
    for i in range(0, len(ids), 900):
        chunk = ids[i:i + 900]
        d = lookup(chunk)
        for k, v in d.items():
            if v:
                res[k] = {'symbol': v.get('display_name'), 'biotype': v.get('biotype'),
                          'description': (v.get('description') or '')[:200]}
        print(f'  batch {i // 900 + 1}: cumulative {len(res)}')
        time.sleep(0.3)
    out = {'source': 'Ensembl REST POST /lookup/id, species bos_taurus, assembly ARS-UCD2.0, '
                     'accessed 2026-10-25',
           'n_queried': len(ids), 'n_resolved': len(res), 'mapping': res}
    (OUT / 'ensembl_symbol_resolution.json').write_text(json.dumps(out, indent=1, ensure_ascii=False),
                                                        encoding='utf-8')
    hi = {r['gene'] for r in rows if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05}
    print(f'resolved {len(res)}/{len(ids)} ({100*len(res)/len(ids):.1f}%)')
    print(f'high-confidence tier resolved: {sum(1 for g in hi if g in res)}/{len(hi)}')


if __name__ == '__main__':
    main()
