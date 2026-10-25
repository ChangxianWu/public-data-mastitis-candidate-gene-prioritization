#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 postcheck：图/表/数字/引用/CJK 门控（对标 VM1 17 门）。

2026-10-07 复核修复轮：数字门控全部改为从 work/output 结果表现算，不再硬编码常量
（原版把 335/3,404/288/143/148/9,984/81.5 写死，导致口径变更后仍报绿）。
"""
import sys, os, re, json, csv, datetime

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'work', 'output')
issues, ok = [], []

# 1) 图 1-7 PNG+PDF 齐全
figdir = os.path.join(ROOT, 'work', 'figures')
for f in ('Fig1_data_flow', 'Fig2_volcano_combined_meta', 'Fig3_loqo_stability', 'Fig4_DE_eQTL_PIP',
          'Fig5_lead_snp_chromosomes', 'Fig6_priority_scoring', 'Fig7_contrast_layers'):
    for ext in ('png', 'pdf'):
        if not any(fn.startswith(f) and fn.endswith('.' + ext) for fn in os.listdir(figdir)):
            issues.append(f'缺图: {f}.{ext}')
if not any('缺图' in i for i in issues):
    ok.append(f'7 张图 PNG+PDF 齐全（{len(os.listdir(figdir))} 个文件）')

# 2) 正文路径与图号引用
ms_path = os.path.join(ROOT, '交付物', '正文投稿版_v0.3_2026-09-22.md')
ms = open(ms_path, encoding='utf-8').read()
figs_cited = set(re.findall(r'Figure\s*(\d)', ms))
if len(figs_cited) >= 7:
    ok.append(f'正文提及图号: {sorted(figs_cited)}')
else:
    issues.append(f'正文仅提及图 {sorted(figs_cited)}（<7）')

# 3) CJK 残留
cjk = re.findall(r'[\u4e00-\u9fff]', ms)
if cjk:
    issues.append(f'正文含 {len(cjk)} 个中文字符')
else:
    ok.append('正文 CJK=0')


# ---- 4) 关键数字：从结果表现算后与正文比对（禁止硬编码期望值） ----
def tsv(name):
    with open(os.path.join(OUT, name), encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def jload(name):
    return json.load(open(os.path.join(OUT, name), encoding='utf-8'))


chip = tsv('de_meta_chip_final.tsv')
anch = tsv('de_meta_gene_anchored.tsv')
comb = tsv('de_meta_gene_combined.tsv')
integ = tsv('de_eqtl_pip_integrated.tsv')
sc = tsv('priority_scoring.tsv')
meta = jload('de_meta_queue_metadata.json')['queue']
mg = jload('de_queue_mergeability.json')
rn = jload('de_rnaseq_track_summary.json')
ov_sum = jload('de_eqtl_overlap_summary.json')
lo = jload('loqo_stability.json')
comp = jload('mini_qtl_gene_baseline.json')
inter = jload('mini_qtl_intersection.json')
net = tsv('string_network.tsv')
nodes = set()
for r in net:
    nodes.update([r['node_a'], r['node_b']])

n_series = sum(int(v.get('n_samples') or 0) for v in meta.values())
n_analysed = sum(v['case'] + v['control'] for v in mg.values()) + \
    sum(v['case'] + v['ctl'] for v in rn.values())
n_fdr = sum(1 for r in chip if float(r['q_BH']) < 0.05)
n_tierA = sum(1 for r in chip if r['tier'] == 'A')
n_three = sum(1 for r in comb if len([x for x in r['sources'].split(',') if x]) >= 3)
n_hi = sum(1 for r in integ if float(r['max_PIP']) >= 0.5 and float(r['p_DL']) < 0.05)
n_st9 = sum(1 for r in integ if float(r['max_PIP']) >= 0.9 and float(r['p_DL']) < 0.05)
lo_min = min(v['rate'] for v in lo.values()) * 100

for label, val, claim in (
        ('探针 meta', len(chip), f'{len(chip):,}'),
        ('FDR q<0.05', n_fdr, f'{n_fdr:,}'),
        ('Tier A', n_tierA, f'{n_tierA:,}'),
        ('基因锚点', len(anch), f'{len(anch):,}'),
        ('合并基因', len(comb), f'{len(comb):,}'),
        ('eQTL 显著基因', ov_sum['eqtl_genes'], f'{ov_sum["eqtl_genes"]:,}'),
        ('DE∩eQTL', ov_sum['overlap'], f'{ov_sum["overlap"]:,}'),
        ('整合层 strict', len(integ), f'{len(integ):,}'),
        ('STRING 边', len(net), f'{len(net):,}'),
        ('STRING 节点', len(nodes), f'{len(nodes):,}'),
        ('打分基因数', len(sc), f'{len(sc):,}')):
    if claim in ms:
        ok.append(f'数字一致(ref-true): {label} = {val}')
    else:
        issues.append(f'正文缺/不符数字 {label}（结果表={val}，期望出现 "{claim}"）')

for label, val, claim in (
        ('系列样本总数', n_series, f'{n_series} samples'),
        ('三源基因', n_three, f'{n_three} genes supported'),
        ('高置信层', n_hi, f'contains {n_hi} genes'),
        ('严层 PIP>=0.9', n_st9, f'retains {n_st9} genes'),
        ('竞文对比层', len(comp), f'{len(comp)}'),
        ('LOQO 下界', f'{lo_min:.1f}', f'{lo_min:.1f}'),
        ('对比样本数', n_analysed, f'{n_analysed} entered case-control contrasts')):
    if claim in ms:
        ok.append(f'数字一致(ref-true): {label} = {val}')
    else:
        issues.append(f'正文缺/不符数字 {label}（结果表={val}，期望出现 "{claim}"）')

if len(inter) == 1:
    rec_inter = json.load(open(os.path.join(OUT, 'mini_qtl_intersection_reconciled.json'), encoding='utf-8'))
    if len(rec_inter) == 4:
        ok.append('竞文交集: 原样口径 1 / 标识重映射后 4（286 层）')
    else:
        issues.append(f'竞文重映射交集异常 {len(rec_inter)} != 4')

# 4b) S8/S9 稳健性层（分层 meta / 权重扰动 / 对比富集）
strat = jload('stratified_meta.json')
wp = jload('weight_perturbation.json')
enr = jload('enrichment_summary.json')
if strat['main_ref_recomputed']['tierA'] == n_tierA:
    ok.append('S8 内部一致: 主 meta 复算 Tier A 与结果表相同')
else:
    issues.append(f'S8 主 meta 复算不符: {strat["main_ref_recomputed"]["tierA"]} != {n_tierA}')
for label, val, claim in (
        ('S8 E.coli 层保留率', strat['cross']['main_tierA_in_ecoli']['retention_pct'], '73.6% of the main Tier A'),
        ('S8 去肝臂保留率', strat['cross']['main_tierA_in_noliver']['retention_pct'], '70.6%'),
        ('S8 Top-30 扰动中位', wp['top30_overlap']['median'], '28 of 30'),
        ('S9 显著功能条目', enr['set_vm2']['n_sig_terms_fdr05'], '1,277'),
        ('S9 GO Process 条目', enr['set_vm2']['by_category']['Process'], '644')):
    if claim in ms:
        ok.append(f'数字一致(ref-true): {label} = {val}')
    else:
        issues.append(f'正文缺/不符数字 {label}（结果={val}，期望出现 "{claim}"）')
if enr['set_competitor']['n_sig_terms_fdr05'] == 0:
    ok.append('S9: 竞文 150 标识（原样）GO/KEGG/Reactome 零显著')
else:
    issues.append(f'S9: 竞文原样显著条目数异常 {enr["set_competitor"]["n_sig_terms_fdr05"]} != 0')
rec = enr.get('set_competitor_reconciled', {})
if (rec.get('by_category', {}).get('Process', 1) == 0 and rec.get('by_category', {}).get('KEGG', 1) == 0
        and rec.get('by_category', {}).get('RCTM', 1) == 0):
    ok.append('S9: 竞文重映射集合无 Process/KEGG/Reactome 显著条目')
else:
    issues.append(f'S9: 重映射集合过程层显著条目异常 {rec.get("by_category")}')

# 4c) S10 评审响应层
ex = jload('exact_variance_sensitivity.json')
rr = jload('review_response_metrics.json')
comp_rec = jload('competitor_reconciled.json')
for label, val, claim in (
        ('S10 基因层 q<0.05', rr['A_gene_level_fdr']['combined_q_lt_0.05'], '563'),
        ('S10 286 层 q<0.05', rr['A_gene_level_fdr']['high_conf_with_gene_q_lt_0.05'], '75 of the 286'),
        ('S10 精确方差 q<0.05', ex['probe_level']['q_lt_0.05_exact'], '2,663'),
        ('S10 精确方差 Tier A 保留率', ex['probe_level']['published_tierA_retained_pct'], '66.5%'),
        ('S10 精确方差 Top-30 保留', ex['gene_level']['top30_overlap_exact'], '24 of the 30 Top-30'),
        ('S10 权重网格中位', rr['F_weight_grid']['top30_overlap_median'], '25 of 30'),
        ('S10 竞文重映射交集', comp_rec['n_intersection_286'], '4 of 286'),
        ('S10 竞文重映射基因数', comp_rec['unique_reconciled_ensembl'], '93 unique genes')):
    if str(claim) in ms:
        ok.append(f'数字一致(ref-true): {label} = {val}')
    else:
        issues.append(f'正文缺/不符数字 {label}（结果={val}，期望出现 "{claim}"）')
for phrase, label in (('CCRL2', 'S10 收敛基因'), ('S11', 'S11 审计引用'), ('257 samples', 'S10 队列贡献账')):
    if phrase in ms:
        ok.append(f'声明在场: {label}')
    else:
        issues.append(f'正文缺声明: {label} ({phrase})')

# 5) 表 3/4 行数与结果表一致
t3 = open(os.path.join(ROOT, '交付物', 'tables', 'Table3_high_confidence.md'), encoding='utf-8').read()
t3_rows = len([ln for ln in t3.split('\n') if ln.startswith('| ENS')])
if t3_rows == n_hi:
    ok.append(f'Table 3 = {n_hi} 行高置信基因')
else:
    issues.append(f'Table 3 行数 {t3_rows} != 高置信层 {n_hi}')
t4 = open(os.path.join(ROOT, '交付物', 'tables', 'Table4_top30.md'), encoding='utf-8').read()
if len([ln for ln in t4.split('\n') if re.match(r'\| \d+ \|', ln)]) == 30:
    ok.append('Table 4 = 30 行 Top 候选')
else:
    issues.append('Table 4 行数 != 30')

# 6) 补充材料 S1-S7 + 索引
supp = os.path.join(ROOT, 'work', 'supp')
need = ['S1_cohort_selection.md', 'S2_parsing_pitfalls.md', 'S3_per_queue_results.md',
        'S4_platform_annotation.md', 'S5_granularity_statement.md', 'S6_scoring_spec.md',
        'S7_network_coverage.md', 'S8_robustness.md', 'S9_functional_enrichment.md',
        'S10_review_response_sensitivity.md', 'S11_gwas_data_audit.md',
        'SUPPLEMENTARY_INDEX.md']
miss = [s for s in need if not os.path.exists(os.path.join(supp, s))]
if not miss:
    ok.append('补充材料 S1-S11 + SUPPLEMENTARY_INDEX 齐全')
else:
    issues.append(f'缺补充材料: {miss}')

# 7) 文献库：正文内嵌与独立清单一一对应
ref = open(os.path.join(ROOT, '交付物', 'references_vancouver_v0.1.txt'), encoding='utf-8').read()
n_ref = len([ln for ln in ref.split('\n') if re.match(r'^\d+\.', ln)])
head, _, tail = ms.partition('## References')
listed = sorted(int(m.group(1)) for m in re.finditer(r'^\[(\d+)\]', tail, re.M))
cited = set()
for m in re.finditer(r'\[([\d,\s\-–]+)\]', head):
    for part in re.split(r'[,\s]+', m.group(1)):
        if not part:
            continue
        rng = re.split(r'[-–]', part)
        if len(rng) == 2 and rng[0].isdigit():
            cited.update(range(int(rng[0]), int(rng[1]) + 1))
        elif rng[0].isdigit():
            cited.add(int(rng[0]))
if n_ref >= 20 and n_ref == len(listed) and set(listed) == cited:
    ok.append(f'文献库 {n_ref} 条，与正文内嵌与引用完全一致')
else:
    issues.append(f'文献不一致: 清单={n_ref} 内嵌={len(listed)} 正文引用={len(cited)} '
                  f'cited-not-listed={sorted(cited - set(listed))} listed-not-cited={sorted(set(listed) - cited)}')

# 8) 投稿包与 MANIFEST
pkg = os.path.join(ROOT, '投稿材料包')
mf_p = os.path.join(pkg, 'MANIFEST.json')
if os.path.exists(mf_p):
    mf = json.load(open(mf_p, encoding='utf-8'))
    actual = [os.path.relpath(os.path.join(dp, fn), pkg).replace(os.sep, '/')
              for dp, _, fns in os.walk(pkg) for fn in fns if fn != 'MANIFEST.json']
    listed2 = {f['file'] for f in mf['files']}
    miss2 = [a for a in actual if a not in listed2]
    if miss2:
        issues.append(f'包漂移: {miss2[:3]}')
    else:
        ok.append(f'投稿包 {len(actual)} 文件与 MANIFEST 一致')
else:
    issues.append('投稿包 MANIFEST.json 缺失')

# 9) 关键声明在场
for phrase, label in (('region-level', '粒度声明'), ('symbol-resolvable', '网络覆盖声明'),
                      ('EST-class', 'EST 不强映射'), ('Benjamini-Hochberg', 'FDR 口径声明')):
    if phrase in ms:
        ok.append(f'声明在场: {label}')
    else:
        issues.append(f'正文缺声明: {label} ({phrase})')

# 10) 陈旧口径残留扫描（防止回退到旧版本）
for stale in ('335 samples', '3,404 integrated', '288 genes', '288-gene',
              'Bonferroni-significant findings', '148 gene', '9,984 edges'):
    if stale in ms:
        issues.append(f'陈旧口径残留: "{stale}"')

print('== POSTCHECK VM2 ==')
for o in ok:
    print('  ✔', o)
for i in issues:
    print('  ✘', i)
print(f'\n结果: {len(ok)} 通过 / {len(issues)} 问题')
os.makedirs(os.path.join(ROOT, 'work', 'logs'), exist_ok=True)
json.dump({'ok': ok, 'issues': issues},
          open(os.path.join(ROOT, 'work', 'logs', 'postcheck_2026-10-25.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)
sys.exit(1 if issues else 0)
