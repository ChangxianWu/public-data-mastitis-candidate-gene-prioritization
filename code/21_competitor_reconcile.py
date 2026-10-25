#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 6c: 竞文标识的权威重映射 + 重算交集 + 重富集（评审 M5）。

15 号脚本按"标识原样"求交：150 个竞文标识中 69 个是 ENSBTAG、81 个是符号；符号
无法直接与 Ensembl 整合层求交，故旧口径只得 1 个共享基因。本脚本把 81 个符号经
Ensembl REST（bos_taurus）映射为 ENSBTAG，与 69 个既有 ENSBTAG 合并去重，再与
286 / 3,395 层重算交集，并以该"已解析集合"重做 STRING 富集（S9 的公平版本）。

输出：
  work/output/competitor_reconciled.json
  work/output/mini_qtl_intersection_reconciled.json
  work/output/enrichment_competitor_reconciled.tsv
  work/output/enrichment_summary.json（追加 set_competitor_reconciled）
运行：python work/code/21_competitor_reconcile.py
"""
import csv, json, os, ssl, sys, time, urllib.parse, urllib.request
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.normpath(os.path.join(HERE, 'output'))
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}), urllib.request.HTTPSHandler(context=CTX))
OPENER.addheaders = [('User-Agent', 'VM2-review/1.0')]


def outpath(name):
    base = os.path.normpath(OUT)
    p = os.path.normpath(base + os.sep + name)
    if not p.startswith(base + os.sep):
        raise SystemExit('output path escaped OUT dir: ' + name)
    return p


def dump(obj, name):
    Path(outpath(name)).write_text(json.dumps(obj, indent=1, ensure_ascii=False), encoding='utf-8')


def rd(name):
    return list(csv.DictReader(open(outpath(name), encoding='utf-8'), delimiter='\t'))


def lookup(url, timeout=60):
    try:
        return json.loads(OPENER.open(url, timeout=timeout).read())
    except Exception:
        return None


def ensembl_symbol(sym):
    d = lookup('https://rest.ensembl.org/lookup/symbol/bos_taurus/'
               + urllib.parse.quote(sym) + '?content-type=application/json')
    return (d or {}).get('id')


def ensembl_name(gid):
    d = lookup('https://rest.ensembl.org/lookup/id/' + urllib.parse.quote(gid)
               + '?content-type=application/json')
    return (d or {}).get('display_name')


def string_enrichment(ids):
    params = urllib.parse.urlencode({'identifiers': '\r'.join(ids), 'species': 9913,
                                     'caller_identity': 'VM2_mastitis'}).encode('utf-8')
    for a in range(4):
        try:
            req = urllib.request.Request('https://string-db.org/api/tsv/enrichment', data=params,
                                         headers={'Content-Type': 'application/x-www-form-urlencoded'})
            body = OPENER.open(req, timeout=180).read().decode('utf-8', 'replace')
            lines = [l for l in body.strip().split('\n') if l.strip()]
            hdr = lines[0].split('\t')
            return [dict(zip(hdr, l.split('\t'))) for l in lines[1:]]
        except Exception as e:
            print('  enrichment retry', a + 1, type(e).__name__)
            time.sleep(2 * (a + 1))
    raise SystemExit('enrichment failed')


def main():
    comp = json.load(open(outpath('mini_qtl_gene_baseline.json'), encoding='utf-8'))
    given_ens = sorted({c for c in comp if c.startswith('ENSBTAG')})
    syms = sorted({c for c in comp if not c.startswith('ENSBTAG')})
    mapping_path = Path(outpath('competitor_symbol_mapping.json'))
    cache = json.loads(mapping_path.read_text(encoding='utf-8')) if mapping_path.is_file() else {}
    sym_map = cache.get('symbol_map', {})
    unresolved = list(cache.get('unresolved', []))
    for s in syms:
        if s in sym_map:
            continue
        gid = ensembl_symbol(s)
        if gid:
            sym_map[s] = {'ensembl': gid, 'assembly': 'ARS-UCD2.0'}
        else:
            if s not in unresolved:
                unresolved.append(s)
        time.sleep(0.05)
    ens2sym = {given: (sym_map.get(given, {}) or {}).get('ensembl') for given in []}
    # 既有 ENSBTAG 取符号（用于富集）
    id2sym = {}
    for gid in given_ens:
        nm = ensembl_name(gid)
        if nm:
            id2sym[gid] = nm
        time.sleep(0.05)
    for s, d in sym_map.items():
        id2sym.setdefault(d['ensembl'], s)
    reconciled = sorted(set(given_ens) | {d['ensembl'] for d in sym_map.values() if d.get('ensembl')})
    hi = [r['gene'] for r in rd('de_eqtl_pip_integrated.tsv')
          if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05]
    integ = [r['gene'] for r in rd('de_eqtl_pip_integrated.tsv')]
    inter_hi = sorted(set(reconciled) & set(hi))
    inter_int = sorted(set(reconciled) & set(integ))
    out = {
        'identifiers_as_published': len(comp),
        'given_ensembl': len(given_ens),
        'given_symbols': len(syms),
        'symbols_resolved': sum(1 for s in syms if s in sym_map),
        'unresolved_symbols': unresolved,
        'unique_reconciled_ensembl': len(reconciled),
        'reconciled_ensembl_with_symbol': len([g for g in reconciled if id2sym.get(g)]),
        'intersection_286_reconciled': inter_hi,
        'intersection_3395_reconciled': inter_int,
        'n_intersection_286': len(inter_hi),
        'n_intersection_3395': len(inter_int),
        'shared_symbols_286': {g: id2sym.get(g) for g in inter_hi},
        'source': 'Ensembl REST bos_taurus lookups (2026-10-25); STRING enrichment on reconciled symbols',
    }
    dump(out, 'competitor_reconciled.json')
    dump(inter_hi, 'mini_qtl_intersection_reconciled.json')
    print(json.dumps(out, indent=1, ensure_ascii=False))

    symbols = sorted({id2sym[g] for g in reconciled if id2sym.get(g)})
    print('enrichment set (reconciled symbols):', len(symbols))
    rows = string_enrichment(symbols)
    cols = ['term', 'category', 'number_of_genes', 'fdr', 'description']
    Path(outpath('enrichment_competitor_reconciled.tsv')).write_text(
        '\t'.join(cols) + '\n' + '\n'.join('\t'.join(str(r.get(c) or '') for c in cols) for r in rows) + '\n',
        encoding='utf-8')
    keep = {'Process', 'Function', 'Component', 'KEGG', 'RCTM'}
    sig = [r for r in rows if r.get('category') in keep and r.get('fdr') and float(r['fdr']) < 0.05]
    by = {c: sum(1 for r in sig if r.get('category') == c) for c in keep}
    summary = json.load(open(outpath('enrichment_summary.json'), encoding='utf-8'))
    summary['set_competitor_reconciled'] = {
        'n_genes': len(symbols), 'n_terms_total': len(rows),
        'n_sig_terms_fdr05': len(sig), 'by_category': by,
        'top_process': [{'term': r.get('description'), 'fdr': r.get('fdr'), 'genes': r.get('number_of_genes')}
                        for r in sorted([x for x in sig if x.get('category') == 'Process'],
                                        key=lambda x: float(x['fdr']))[:8]],
    }
    dump(summary, 'enrichment_summary.json')
    print(json.dumps(summary['set_competitor_reconciled'], indent=1, ensure_ascii=False)[:1500])


if __name__ == '__main__':
    main()
