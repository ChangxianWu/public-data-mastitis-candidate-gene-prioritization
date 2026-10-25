#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 4d: 功能富集对比层（S9 补充材料的数据源）。

  集合 A：本研究 764 三源验证基因（chip anchors + 双 RNA-seq 队列，sources>=3；全部有 symbol）
  集合 B：竞文 PMID 39009986 的 150 个 gene-based 命中（mini_qtl_gene_baseline.json，symbol）
  端点：STRING v11.5 /api/tsv/enrichment（species 9913，全基因组背景，一次性整集提交）
  目的：两条证据流是否收敛于同一生物学过程（贡献"收敛-分歧"解释层）。

输出 work/output/enrichment_vm2_764.tsv + enrichment_competitor_150.tsv +
enrichment_summary.json。证书过期 → unverified SSL 上下文（与 09 号一致）。
"""
import csv, json, os, ssl, sys, time, urllib.parse, urllib.request
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.normpath(os.path.join(HERE, 'output'))


def outpath(name):
    """OUT 目录下的静态文件名拼接：归一化并校验不越出 OUT（安全门要求）。"""
    base = os.path.normpath(OUT)
    p = os.path.normpath(base + os.sep + name)
    if not p.startswith(base + os.sep):
        raise SystemExit('output path escaped OUT dir: ' + name)
    return p


def dump_json(obj, name):
    Path(outpath(name)).write_text(json.dumps(obj, indent=1), encoding='utf-8')


CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}), urllib.request.HTTPSHandler(context=CTX))
OPENER.addheaders = [('User-Agent', 'VM2-gate/1.0')]

KEEP_CATEGORIES = {'Process', 'Function', 'Component', 'KEGG', 'RCTM'}


def term_key(r):
    """term_id 列在 v11.5 返回中常为空；以 term（如 GOBP:0006955）作条目身份。"""
    return r.get('term') or r.get('term_id')


def string_enrichment(ids, tag):
    """整集一次性提交（富集是集合级统计，禁止分批）；4 次重试 + 退避。"""
    params = urllib.parse.urlencode({'identifiers': '\r'.join(ids), 'species': 9913,
                                     'caller_identity': 'VM2_mastitis'}).encode('utf-8')
    last = None
    for attempt in range(4):
        try:
            req = urllib.request.Request('https://string-db.org/api/tsv/enrichment', data=params,
                                         headers={'User-Agent': 'VM2-gate/1.0',
                                                  'Content-Type': 'application/x-www-form-urlencoded'})
            body = OPENER.open(req, timeout=180).read().decode('utf-8', 'replace')
            lines = [l for l in body.strip().split('\n') if l.strip()]
            hdr = lines[0].split('\t')
            rows = [dict(zip(hdr, l.split('\t'))) for l in lines[1:]]
            print(f'  {tag}: {len(ids)} identifiers -> {len(rows)} enrichment rows')
            return rows
        except Exception as e:  # 端点抖动重试
            last = e
            print(f'  {tag}: attempt {attempt + 1} failed: {type(e).__name__}: {str(e)[:80]}')
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f'enrichment failed for {tag}: {last}')


def write_tsv(rows, name):
    cols = ['term_id', 'term', 'category', 'number_of_genes', 'fdr', 'description', 'preferredNames']
    lines = ['\t'.join(cols)]
    for r in rows:
        lines.append('\t'.join(str(r.get(k) or '') for k in cols))
    Path(outpath(name)).write_text('\n'.join(lines) + '\n', encoding='utf-8')


def top_terms(rows, cat, n):
    sel = [r for r in rows if r.get('category') == cat and float(r.get('fdr') or 1) < 0.05]
    sel.sort(key=lambda r: float(r['fdr']))
    return [{'term': r.get('description') or r.get('term'), 'fdr': r.get('fdr'),
             'genes': r.get('number_of_genes')} for r in sel[:n]]


def main():
    comb = list(csv.DictReader(open(outpath('de_meta_gene_combined.tsv'), encoding='utf-8'), delimiter='\t'))
    three = [r['symbol'] for r in comb if len([x for x in r['sources'].split(',') if x]) >= 3 and r.get('symbol')]
    three = sorted(set(three))
    comp = sorted(json.load(open(outpath('mini_qtl_gene_baseline.json'), encoding='utf-8')))

    print('STRING enrichment (one-shot set-level calls):')
    rows_a = string_enrichment(three, 'VM2-764')
    rows_b = string_enrichment(comp, 'competitor-150')

    write_tsv(rows_a, 'enrichment_vm2_764.tsv')
    write_tsv(rows_b, 'enrichment_competitor_150.tsv')

    sig_a = {term_key(r): r for r in rows_a if r.get('category') in KEEP_CATEGORIES and float(r.get('fdr') or 1) < 0.05}
    sig_b = {term_key(r): r for r in rows_b if r.get('category') in KEEP_CATEGORIES and float(r.get('fdr') or 1) < 0.05}
    shared_terms = sorted(set(sig_a) & set(sig_b))

    summary = {
        'endpoint': 'STRING v11.5 /api/tsv/enrichment, species 9913, genome-wide background, one-shot set-level call',
        'set_vm2': {'n_genes': len(three), 'n_sig_terms_fdr05': len(sig_a),
                    'by_category': {c: sum(1 for r in sig_a.values() if r['category'] == c) for c in KEEP_CATEGORIES},
                    'top_process': top_terms(rows_a, 'Process', 8),
                    'top_kegg': top_terms(rows_a, 'KEGG', 8)},
        'set_competitor': {'n_genes': len(comp), 'n_sig_terms_fdr05': len(sig_b),
                           'by_category': {c: sum(1 for r in sig_b.values() if r['category'] == c) for c in KEEP_CATEGORIES},
                           'top_process': top_terms(rows_b, 'Process', 8),
                           'top_kegg': top_terms(rows_b, 'KEGG', 8)},
        'shared_sig_terms': len(shared_terms),
        'shared_term_ids': shared_terms[:20],
        'competitor_go_kegg_sig': {c: sum(1 for r in sig_b.values() if r['category'] == c) for c in ('Process', 'Function', 'Component', 'KEGG', 'RCTM')},
        'note': 'the 286-gene high-confidence tier is deliberately not enriched (10/286 symbol-resolvable, S7); the 764 tri-validated layer is the symbol-complete interpretability surface',
    }
    dump_json(summary, 'enrichment_summary.json')
    print(json.dumps(summary, indent=1, ensure_ascii=False)[:2600])


if __name__ == '__main__':
    main()
