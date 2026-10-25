#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""External replication of the mastitis expression signature in held-out cohorts.

Two cohorts that were screened but not used in the meta-analysis:

  A. GSE94056 (GPL17645 array; 12 cows; E. coli or S. aureus infected quarters at 1-3 h vs an
     uninfected control quarter; 49 arrays) - early in vivo response.
  B. GSE334132 (RNA-seq, MAC-T bovine mammary epithelial cells; parental S. aureus and its small
     colony variant vs uninfected wild-type cells; 15 libraries, TPM) - in vitro, all genes.

Tests, for genes measurable in each external cohort:
  1. direction concordance of the 286 high-confidence screening genes, the 75-gene FDR core and
     the Top-30 with the external infection effect (binomial test against 0.5);
  2. Fisher enrichment of those sets among the external differentially expressed genes;
  3. the five canonical axis genes.

Output: work/output/external_validation.json
Run: python work/code/24_external_validation.py
"""
import csv, gzip, json, re, sys, tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = Path(__file__).resolve().parents[1]
OUT = HERE / 'output'
EXT = HERE / 'data' / 'external'
BRIDGE = HERE / 'data' / 'platform' / 'entrez2ensembl_bridge.json'


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    q = [1.0] * n
    run = 1.0
    for rank, i in enumerate(reversed(order), 1):
        run = min(run, ps[i] * n / max(n - rank + 1, 1))
        q[i] = run
    return q


def welch(mat, case_idx, ctl_idx):
    x1, x0 = mat[:, case_idx], mat[:, ctl_idx]
    n1, n0 = x1.shape[1], x0.shape[1]
    d = np.nanmean(x1, axis=1) - np.nanmean(x0, axis=1)
    se = np.sqrt(np.nanvar(x1, axis=1, ddof=1) / n1 + np.nanvar(x0, axis=1, ddof=1) / n0)
    se[se == 0] = np.nan
    t = d / se
    p = 2 * stats.t.sf(np.abs(t), n1 + n0 - 2)
    return d, p


# ----------------------------------------------------------------------------- cohort A
def cohort_gse94056():
    txt = gzip.open(EXT / 'GSE94056_series_matrix.txt.gz', 'rt', encoding='utf-8', errors='replace').read()
    lines = txt.splitlines()
    titles = [v.strip('"') for v in next(ln for ln in lines if ln.startswith('!Sample_title')).split('\t')[1:]]
    i = next(k for k, ln in enumerate(lines) if 'series_matrix_table_begin' in ln)
    probes, rows = [], []
    for ln in lines[i + 2:]:
        if 'series_matrix_table_end' in ln:
            break
        parts = ln.split('\t')
        probes.append(parts[0].strip('"'))
        rows.append([float(x) if x.strip('"') not in ('', 'null') else np.nan for x in parts[1:]])
    mat = np.array(rows, dtype=float)

    groups = []
    for t in titles:
        tl = t.lower()
        if 'uninfected' in tl or 'control' in tl:
            groups.append('control')
        elif re.search(r's\.\s*aureus', tl):
            groups.append('saur')
        elif re.search(r'e\.\s*coli', tl):
            groups.append('ecoli')
        else:
            groups.append('other')

    p2e = {}
    s = gzip.open(EXT / 'GPL17645_family.soft.gz', 'rt', encoding='utf-8', errors='replace').read()
    started = False
    for ln in s.splitlines():
        if 'platform_table_begin' in ln:
            started = True; continue
        if not started:
            continue
        if 'platform_table_end' in ln:
            break
        parts = ln.split('\t')
        if len(parts) >= 2:
            m = re.search(r'GeneID:(\d+)', parts[1])
            if m:
                p2e[parts[0]] = m.group(1)
    bridge = json.load(open(BRIDGE, encoding='utf-8'))

    by_gene = defaultdict(list)
    for k, pr in enumerate(probes):
        e = p2e.get(pr)
        if e and e in bridge:
            by_gene[bridge[e]].append(k)
    genes, rows = [], []
    for ens, idxs in by_gene.items():
        sub = mat[idxs, :]
        best = idxs[int(np.nanargmax(np.nanmean(sub, axis=1)))]
        genes.append(ens); rows.append(mat[best, :])
    x = np.array(rows)

    ctl = [k for k, g in enumerate(groups) if g == 'control']
    ec = [k for k, g in enumerate(groups) if g == 'ecoli']
    sa = [k for k, g in enumerate(groups) if g == 'saur']
    contrasts = {}
    for lab, case in (('all', ec + sa), ('ecoli', ec), ('saur', sa)):
        d, p = welch(x, case, ctl)
        contrasts[lab] = {'d': d, 'p': p,
                          'q': np.array(bh([min(v, 1.0) if np.isfinite(v) else 1.0 for v in p]))}
    return {'id': 'GSE94056',
            'design': '12 cows; E. coli or S. aureus quarters at 1-3 h vs uninfected quarter; GPL17645 array',
            'n_arrays': int(mat.shape[1]), 'groups': {k: groups.count(k) for k in set(groups)},
            'genes': genes, 'contrasts': contrasts}


# ----------------------------------------------------------------------------- cohort B
def cohort_gse334132():
    import openpyxl
    title_groups = {'MAC-T': 'control', 'MACT-PS': 'saur', 'MACT-SCV': 'saur',
                    'TLR4KO-MACT-PS': 'skip', 'TLR4KO-MACT-SCV': 'skip'}
    data = {}
    with tarfile.open(EXT / 'GSE334132_RAW.tar') as tf:
        for member in tf.getnames():
            if not member.endswith('.xlsx'):
                continue
            tmp = EXT / '_tmp.xlsx'
            tmp.write_bytes(tf.extractfile(member).read())
            wb = openpyxl.load_workbook(tmp, read_only=True)
            ws = wb.active
            for k, row in enumerate(ws.iter_rows(values_only=True)):
                if k == 0:
                    continue
                gid = row[3] if len(row) > 3 else None
                tpm = row[1] if len(row) > 1 else None
                if gid:
                    data.setdefault(str(gid), {})[member] = float(tpm or 0)
            wb.close()
            tmp.unlink(missing_ok=True)

    txt = gzip.open(EXT / 'GSE334132_series_matrix.txt.gz', 'rt', encoding='utf-8', errors='replace').read()
    gsms = [v.strip('"') for v in re.search(r'^!Sample_geo_accession\s+(.*)$', txt, re.M).group(1).split('\t')]
    titles = [v.strip('"') for v in re.search(r'^!Sample_title\s+(.*)$', txt, re.M).group(1).split('\t')]
    gsm2group = {gsm: title_groups.get(re.sub(r'-\d+$', '', t), 'skip') for gsm, t in zip(gsms, titles)}
    members = sorted({m for v in data.values() for m in v})
    groups = [gsm2group.get(m.split('_')[0], 'skip') for m in members]
    genes = sorted(data)
    x = np.log2(np.array([[data[g].get(m, np.nan) for m in members] for g in genes], dtype=float) + 1.0)
    ctl = [k for k, g in enumerate(groups) if g == 'control']
    sa = [k for k, g in enumerate(groups) if g == 'saur']
    d, p = welch(x, sa, ctl)
    q = np.array(bh([min(v, 1.0) if np.isfinite(v) else 1.0 for v in p]))
    return {'id': 'GSE334132',
            'design': 'MAC-T cells; S. aureus parental + SCV vs uninfected wild type; RNA-seq TPM',
            'n_arrays': len(members), 'groups': {k: groups.count(k) for k in set(groups)},
            'genes': genes, 'contrasts': {'saur': {'d': d, 'p': p, 'q': q}}}


# ----------------------------------------------------------------------------- cohort C
def cohort_gse133333():
    """M. bovis infected vs non-infected mammary tissue (very small: 2 vs 3)."""
    import openpyxl
    from collections import defaultdict as dd
    wb = openpyxl.load_workbook(EXT / 'GSE133333_mRNA_FPKM.xlsx', read_only=True)
    ws = wb.active
    it = ws.iter_rows(values_only=True)
    hdr = list(next(it))
    cols = {c: i for i, c in enumerate(hdr) if c in ('infect2', 'infect3', 'cont1', 'cont2', 'cont3')}
    order = [c for c in ('infect2', 'infect3', 'cont1', 'cont2', 'cont3') if c in cols]
    data = dd(list)
    for row in it:
        gid = row[1]
        if not gid:
            continue
        data[gid].append([float(row[cols[c]] or 0) for c in order])
    wb.close()
    genes = list(data)
    x = np.log2(np.array([np.mean(v, axis=0) for v in data.values()], dtype=float) + 1.0)
    case = [k for k, c in enumerate(order) if c.startswith('infect')]
    ctl = [k for k, c in enumerate(order) if c.startswith('cont')]
    d, p = welch(x, case, ctl)
    q = np.array(bh([min(v, 1.0) if np.isfinite(v) else 1.0 for v in p]))
    return {'id': 'GSE133333',
            'design': 'M. bovis infected vs non-infected bovine mammary tissue; RNA-seq FPKM (2 vs 3)',
            'n_arrays': len(order), 'groups': {'infected': len(case), 'control': len(ctl)},
            'genes': genes, 'contrasts': {'saur': {'d': d, 'p': p, 'q': q}}}


# ----------------------------------------------------------------------------- tests
def evaluate(cohort, sets, axis_sym2ens):
    genes = cohort['genes']
    gidx = {g: i for i, g in enumerate(genes)}
    bg = set(genes)
    out = {'id': cohort['id'], 'design': cohort['design'], 'n_arrays': cohort['n_arrays'],
           'groups': cohort['groups'], 'n_genes_external': len(genes), 'contrasts': {}}

    def concord(rows, lab):
        ext = cohort['contrasts'][lab]['d']; q = cohort['contrasts'][lab]['q']
        # genes with an exactly zero external effect (e.g. zero TPM in both groups) carry no
        # direction and must not be counted as discordant
        pairs = [(r, gidx[r['gene']]) for r in rows
                 if r['gene'] in gidx and np.isfinite(ext[gidx[r['gene']]]) and ext[gidx[r['gene']]] != 0]
        if not pairs:
            return None
        same = sum(1 for r, i in pairs if np.sign(ext[i]) == np.sign(float(r['_beta'])))
        n = len(pairs)
        sig = [(r, i) for r, i in pairs if q[i] < 0.05]
        same_sig = sum(1 for r, i in sig if np.sign(ext[i]) == np.sign(float(r['_beta'])))
        return {'n_mapped': n, 'same_sign': same, 'pct': round(100 * same / n, 1),
                'binomial_p': float(f'{stats.binomtest(same, n, 0.5).pvalue:.3g}'),
                'n_external_sig': len(sig), 'same_sign_among_external_sig': same_sig,
                'pct_among_sig': round(100 * same_sig / len(sig), 1) if sig else None}

    def fisher(rows, lab):
        sel = {r['gene'] for r in rows if r['gene'] in gidx}
        sig = {g for g, qq in zip(genes, cohort['contrasts'][lab]['q']) if qq < 0.05}
        a = len(sel & sig); b = len(sel) - a; c = len(sig) - a; d = len(bg) - a - b - c
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        return {'n_selected_mapped': len(sel), 'overlap_with_external_sig': a,
                'expected': round(len(sel) * len(sig) / len(bg), 1),
                'odds_ratio': (round(float(odds), 2) if np.isfinite(odds) else None),
                'fisher_p': float(f'{p:.3g}')}

    for lab in cohort['contrasts']:
        out['contrasts'][lab] = {
            'external_q05': int((cohort['contrasts'][lab]['q'] < 0.05).sum()),
            'three_source_764': concord(sets['three'], lab),
            'high_confidence_286': concord(sets['hi'], lab),
            'fdr_core_75': concord(sets['core'], lab),
            'top30': concord(sets['top30'], lab),
            'enrichment_764': fisher(sets['three'], lab),
            'enrichment_286': fisher(sets['hi'], lab),
            'enrichment_core': fisher(sets['core'], lab),
        }
    out['axis_genes'] = {}
    lab = 'all' if 'all' in cohort['contrasts'] else list(cohort['contrasts'])[0]
    for sym, ens_list in axis_sym2ens.items():
        hit = [g for g in ens_list if g in gidx]
        if not hit:
            out['axis_genes'][sym] = {'measured_externally': False}
            continue
        i = gidx[hit[0]]
        out['axis_genes'][sym] = {'measured_externally': True,
                                  'external_logfc': round(float(cohort['contrasts'][lab]['d'][i]), 3),
                                  'external_q': float(f'{cohort["contrasts"][lab]["q"][i]:.3g}')}
    return out


def main():
    integ = list(csv.DictReader(open(OUT / 'de_eqtl_pip_integrated.tsv', encoding='utf-8'), delimiter='\t'))
    for r in integ:
        r['_beta'] = r['beta_meta']
    hi = [r for r in integ if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05]
    defined = [r for r in integ if r['p_DL'] not in ('', 'NA')]
    qmap = {r['gene']: q for r, q in zip(defined, bh([float(r['p_DL']) for r in defined]))}
    core = [r for r in hi if qmap.get(r['gene'], 1.0) < 0.05]
    sc = list(csv.DictReader(open(OUT / 'priority_scoring.tsv', encoding='utf-8'), delimiter='\t'))[:30]
    for r in sc:
        r['_beta'] = r['beta']
    comb = list(csv.DictReader(open(OUT / 'de_meta_gene_combined.tsv', encoding='utf-8'), delimiter='	'))
    three = [r for r in comb if len([x for x in r['sources'].split(',') if x]) >= 3]
    for r in three:
        r['_beta'] = r['beta_meta']
    sets = {'three': three, 'hi': hi, 'core': core, 'top30': sc}

    bridge = json.load(open(BRIDGE, encoding='utf-8'))
    anchors = list(csv.DictReader(open(OUT / 'de_meta_gene_anchored.tsv', encoding='utf-8'), delimiter='\t'))
    axis = {sym: [bridge[r['entrez_id']] for r in anchors
                  if r['gene_symbol'] == sym and r['entrez_id'] in bridge]
            for sym in ('SAA3', 'CCL2', 'SOCS3', 'TIMP1', 'SERPINE1')}

    cohorts = [cohort_gse94056(), cohort_gse334132(), cohort_gse133333()]
    out = {'layers': {'high_confidence_286': len(hi), 'fdr_core_75': len(core), 'top30': len(sc)},
           'cohorts': [evaluate(c, sets, axis) for c in cohorts]}
    (OUT / 'external_validation.json').write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding='utf-8')
    for c in out['cohorts']:
        print('==', c['id'], '|', c['design'])
        for lab, v in c['contrasts'].items():
            t764 = v['three_source_764']; h = v['high_confidence_286']; k = v['fdr_core_75']; e = v['enrichment_286']
            print(f"  [{lab}] external q<0.05={v['external_q05']} | 764 same={t764['same_sign']}/{t764['n_mapped']} ({t764['pct']}%, p={t764['binomial_p']}) | 286 mapped={h['n_mapped']} "
                  f"same={h['same_sign']} ({h['pct']}%, p={h['binomial_p']}) | core {k['same_sign']}/{k['n_mapped']} | "
                  f"enrich {e['overlap_with_external_sig']}/{e['n_selected_mapped']} (exp {e['expected']}, p={e['fisher_p']})")
        print('  axis:', json.dumps(c.get('axis_genes', {}), ensure_ascii=False))


if __name__ == '__main__':
    main()
