#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 5b: 正文表 1-4 生成器（复核修复轮补齐；此前表 md 无落盘脚本）。

运行方式：在 VM2 课题根目录执行  python work/code/16_make_tables.py

口径（2026-10-07 起 strict）：整合层 = DE∩eQTL显著∩finemapping = 3,395；
高置信层 = PIP>=0.5 & p<0.05 = 286；三源基因 = 764（sources 含 chip+两 RNA-seq）。

输入：work/output/ 各结果表 + de_meta_queue_metadata.json
输出：交付物/tables/Table1_cohorts.md / Table2_meta_summary.md /
      Table3_high_confidence.md（全 286 行）/ Table3_high_confidence_top10.md /
      Table4_top30.md
Table 1 的 Notes 列为编辑性注释（历史定稿文案，随本脚本走，改注释改这里）。
"""
import csv, json, math, os, sys
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
OUTDIR = os.path.abspath('work/output')
TABDIR = os.path.abspath(os.path.join('交付物', 'tables'))
for d in (OUTDIR, TABDIR):
    if not os.path.isdir(d):
        sys.exit(f'目录不存在（请在 VM2 课题根目录运行）: {d}')


def rd(name):
    return os.path.join(OUTDIR, name)


def w(path, text):
    Path(path).write_text(text, encoding='utf-8')
    print('wrote', os.path.basename(path))


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    q = [1.0] * n
    run = 1.0
    for rank, i in enumerate(reversed(order), 1):
        run = min(run, ps[i] * n / max(n - rank + 1, 1))
        q[i] = run
    return q


def main():
    meta = json.load(open(rd('de_meta_queue_metadata.json'), encoding='utf-8'))['queue']
    design = list(csv.DictReader(open(rd('de_queue_design.tsv'), encoding='utf-8'), delimiter='\t'))
    from collections import Counter
    cc = Counter((r['gse'], r['group']) for r in design)
    rn = json.load(open(rd('de_rnaseq_track_summary.json'), encoding='utf-8'))

    # ---- Table 1 ----
    notes = {
        'GSE93082': 'field cases; endometritis excluded',
        'GSE131607': 'milk somatic cells; RPKM->log2',
        'GSE75379': 'udder 24h E.coli vs control; liver track separate',
    }
    contrib = {
        'GSE24560': 'probe meta', 'GSE25319': 'probe meta (liver arm)',
        'GSE24217': 'probe meta (udder arm; same cows as GSE25319)',
        'GSE15025': 'probe meta', 'GSE93082': 'probe meta',
        'GSE131607': 'gene meta (RNA-seq)', 'GSE75379': 'gene meta (RNA-seq; udder arm)',
        'GSE9685': 'per-queue only (matrix lacks probe IDs)',
        'GSE34031': 'per-queue only (matrix lacks probe IDs)',
    }
    lines = ['| Accession | Platform | Series type | Samples | Case | Control | Contribution | Notes |',
             '|---|---|---|---|---|---|---|---|']
    TYPE_LABEL = {'Expression profiling by array': 'Microarray',
                  'Expression profiling by high throughput sequencing': 'RNA-seq'}
    for g in sorted(meta):
        v = meta[g]
        raw_type = v.get('!Series_type', '')
        st = TYPE_LABEL.get(raw_type.strip(), raw_type[:38])
        if g in rn:
            case, ctl = rn[g]['case'], rn[g]['ctl']
        else:
            case, ctl = cc[(g, 'case')], cc[(g, 'control')]
        lines.append(f"| {g} | {v.get('!Series_platform_id','')} | {st} | {v.get('n_samples')} "
                     f"| {case} | {ctl} | {contrib.get(g, '')} | {notes.get(g, '')} |")
    lines.append('')
    lines.append('Nine cohorts (339 samples) were assembled; quantitative synthesis used seven '
                 '(257 samples: five GPL2112 microarray cohorts contributing 233 case-control '
                 'samples to the probe-level meta, and two RNA-seq cohorts contributing 24). '
                 'GSE9685 and GSE34031 were analysed per queue but cannot enter the probe-level '
                 'meta because their deposited series matrices carry numeric row indices rather '
                 'than probe identifiers. GSE75379 contributes its udder arm only (24 liver '
                 'samples excluded); six endometritis samples of GSE93082 are excluded as a '
                 'non-mastitis phenotype.')
    w(os.path.join(TABDIR, 'Table1_cohorts.md'), '\n'.join(lines) + '\n')

    # ---- Table 2 ----
    chip = list(csv.DictReader(open(rd('de_meta_chip_final.tsv'), encoding='utf-8'), delimiter='\t'))
    q05 = sum(1 for r in chip if float(r['q_BH']) < 0.05)
    ta = sum(1 for r in chip if r['tier'] == 'A')
    tb = sum(1 for r in chip if r['tier'] == 'B')
    tc = sum(1 for r in chip if r['tier'] == 'C')
    anch = list(csv.DictReader(open(rd('de_meta_gene_anchored.tsv'), encoding='utf-8'), delimiter='\t'))
    up = sum(1 for r in anch if float(r['beta_meta']) > 0)
    dn = len(anch) - up
    comb = list(csv.DictReader(open(rd('de_meta_gene_combined.tsv'), encoding='utf-8'), delimiter='\t'))
    cp05 = sum(1 for r in comb if r['p_DL'] and float(r['p_DL']) < 0.05)
    three = sum(1 for r in comb if len([x for x in r['sources'].split(',') if x]) >= 3)
    ov = [l.split('\t') for l in open(rd('de_eqtl_overlap_genes.txt'), encoding='utf-8')]
    ov05 = sum(1 for c in ov if c[3] and float(c[3]) < 0.05)
    integ = list(csv.DictReader(open(rd('de_eqtl_pip_integrated.tsv'), encoding='utf-8'), delimiter='\t'))
    hi = sum(1 for r in integ if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05)
    # 基因层 BH（评审 M1）：combined 以有定义 p 的基因为检验集；integrated 以整合层为检验集
    comb_def = [(r, float(r['p_DL'])) for r in comb if r['p_DL'] and math.isfinite(float(r['p_DL']))]
    comb_q = bh([p for _, p in comb_def])
    comb_q05 = sum(1 for q in comb_q if q < 0.05)
    comb_p05 = sum(1 for _, p in comb_def if p < 0.05)
    int_def = [(r, float(r['p_DL'])) for r in integ if math.isfinite(float(r['p_DL']))]
    int_q = bh([p for _, p in int_def])
    int_qmap = {r['gene']: q for (r, _), q in zip(int_def, int_q)}
    int_q05 = sum(1 for q in int_q if q < 0.05)
    hi_q05 = sum(1 for r in integ if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05
                 and int_qmap.get(r['gene'], 1.0) < 0.05)
    t2 = ['| Layer | N | BH q<0.05 | Direction-consistent |', '|---|---|---|---|',
          f'| Probe-level (>=3 queues; five GPL2112 cohorts) | {len(chip):,} | {q05:,} | Tier A {ta:,} / B {tb:,} / C {tc:,} |',
          f'| Gene-anchored (symbol layer) | {len(anch):,} | - | up {up} / down {dn} |',
          f'| Combined gene-level | {len(comb):,} | {comb_q05:,} (of {len(comb_def):,} with defined p; {comb_p05:,} at nominal p<0.05) | three-source {three} |',
          f'| x mammary eQTL overlap | {len(ov):,} | - | - |',
          f'| Integrated (DE x significant cis-eQTL x fine-mapped) | {len(integ):,} | {int_q05:,} (of {len(int_def):,}) | high-confidence {hi} (q<0.05: {hi_q05}) |',
          '',
          'Notes. Gene-anchored is a symbol-level layer (chip probes with a platform symbol/Entrez '
          'anchor, one row per gene at its most significant probe); the combined gene-level layer '
          'maps anchors to Ensembl and adds the two RNA-seq cohorts. The two layers are therefore '
          'nested but not identical, and the 764 three-source genes are defined on the combined layer. '
          'Gene-level BH is reported as a sensitivity (M1): the high-confidence tier is a PIP-based '
          'screening tier, not an FDR-controlled tier; all 30 Top-30 candidates pass gene-level '
          'q<0.05. Under the exact Welch-variance convention (S10), probe-level q<0.05 = 2,663 and '
          '66.5% of the published Tier A probes remain significant with identical direction.']
    w(os.path.join(TABDIR, 'Table2_meta_summary.md'), '\n'.join(t2) + '\n')

    # ---- Table 3（全 286 行 + top10）----
    hi_rows = [r for r in integ if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05]
    hi_rows.sort(key=lambda r: float(r['p_DL']))

    def t3_row(r):
        q = int_qmap.get(r['gene'])
        qs = f'{q:.2e}' if q is not None else 'NA'
        return (f"| {r['gene']} | {r['symbol']} | {r['n_queues']} | {float(r['beta_meta']):+.2f} "
                f"| {float(r['p_DL']):.2e} | {qs} | {float(r['max_PIP']):.2f} | {r['lead_SNP']} |")

    t3 = ['| Gene | Symbol | Queues | Beta | p (meta) | q (BH, integrated) | Max PIP | Lead SNP |',
          '|---|---|---|---|---|---|---|---|']
    t3 += [t3_row(r) for r in hi_rows]
    rec = json.load(open(rd('competitor_reconciled.json'), encoding='utf-8'))
    shared_txt = '; '.join(f"{sym} ({gid})" for gid, sym in sorted(rec['shared_symbols_286'].items()))
    t3.append('')
    t3.append('Shared with the competitor gene-based set after identifier reconciliation: '
              + shared_txt + '. ENSBTAG00000006155 (CCRL2) is rank 4 in Table 4.')
    w(os.path.join(TABDIR, 'Table3_high_confidence.md'), '\n'.join(t3) + '\n')
    t3t = ['| Gene | Symbol | Queues | Beta | p (meta) | q (BH, integrated) | Max PIP | Lead SNP |',
           '|---|---|---|---|---|---|---|---|']
    t3t += [t3_row(r) for r in hi_rows[:10]]
    t3t.append('')
    t3t.append('Shared with the competitor gene-based set after identifier reconciliation: '
               + shared_txt + '. ENSBTAG00000006155 (CCRL2) is rank 4 in Table 4.')
    w(os.path.join(TABDIR, 'Table3_high_confidence_top10.md'), '\n'.join(t3t) + '\n')

    # ---- Table 4（Top-30，直接排版打分表）----
    top = list(csv.DictReader(open(rd('top_candidates.tsv'), encoding='utf-8'), delimiter='\t'))[:30]
    t4 = ['| Rank | Gene | Symbol | Score | S1 DE | S3 PIP | S4 Net | STRING degree | Queues | Beta | p | Lead SNP | Validation |',
          '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in top:
        t4.append(f"| {r['rank']} | {r['gene']} | {r['symbol']} | {r['score']} | {round(float(r['S1_DE']),3)} "
                  f"| {round(float(r['S3_PIP']),3)} | {round(float(r['S4_network']),3)} | {r.get('STRING_degree', 0)} "
                  f"| {r['n_queues']} | {float(r['beta']):+.2f} | {float(r['p_DL']):.1e} | {r['lead_SNP']} "
                  f"| {r['validation_suggestion']} |")
    t4.append('')
    t4.append('Score = 0.35 x S1 + 0.25 x S2 + 0.30 x S3 + 0.10 x S4, where S1 is the min-max '
              'normalized -log10 combined meta p, S2 = 0.5 for every integrated gene (eQTL presence '
              'is an inclusion criterion, constant within the tier and therefore rank-inactive), '
              'S3 = maximum fine-mapping PIP, and S4 = normalized STRING subnetwork weighted degree '
              '(zero for genes without a symbol-resolvable STRING hit; 11 of 3,392 scored genes '
              'have a non-zero S4). STRING degree is the unweighted degree in the frozen seed '
              'subnetwork. The network component is an annotation layer, not a ranking driver: '
              'only CA6 among the 30 entries has a non-zero S4 (contribution 0.0035); a ±50% '
              'weight perturbation retains a median of 28 of 30 entries and a 103-point grid '
              'retains a median of 25 of 30 (S8, S10). Validation suggestions are stratified by '
              'executability: genes without an official symbol receive lead-SNP genotyping and '
              'annotation-dependent follow-up, not qPCR primer suggestions. Rank 4 '
              '(ENSBTAG00000006155) is CCRL2, one of the four genes shared with the competitor '
              'after identifier reconciliation (S6).')
    w(os.path.join(TABDIR, 'Table4_top30.md'), '\n'.join(t4) + '\n')

    # ---- Table 4（GSE 投稿版：9 列，去掉长文案列，避免排版挤压）----
    t4m = ['| Rank | Gene | Symbol | Score | S1 DE | S3 PIP | S4 Net | p | Lead SNP |',
           '|---|---|---|---|---|---|---|---|---|']
    for r in top:
        t4m.append(f"| {r['rank']} | {r['gene']} | {r['symbol']} | {r['score']} | {round(float(r['S1_DE']),3)} "
                   f"| {round(float(r['S3_PIP']),3)} | {round(float(r['S4_network']),3)} "
                   f"| {float(r['p_DL']):.1e} | {r['lead_SNP']} |")
    t4m.append('')
    t4m.append('Score = 0.35 x S1 + 0.25 x S2 + 0.30 x S3 + 0.10 x S4, where S1 is the min-max '
               'normalized -log10 combined meta-analysis p, S2 is eQTL presence (0.5 for every '
               'integrated gene and therefore rank-inactive), S3 is the maximum fine-mapping PIP, '
               'and S4 is the normalized STRING subnetwork weighted degree (zero for genes without '
               'a symbol-resolvable STRING hit). Per-gene follow-up suggestions, the STRING '
               'subnetwork degree and the unweighted degree are provided in Additional file 2 '
               '(Table S11). Rank 4 (ENSBTAG00000006155) is CCRL2, one of the four '
               'genes shared with the competitor after identifier reconciliation.')
    w(os.path.join(TABDIR, 'Table4_top30_manuscript.md'), '\n'.join(t4m) + '\n')
    print(f'Table3 rows={len(hi_rows)}  Table4 rows={len(top)}')


if __name__ == '__main__':
    main()
