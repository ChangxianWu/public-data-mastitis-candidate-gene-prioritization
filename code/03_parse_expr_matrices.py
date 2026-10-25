#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 1b: 表达矩阵解析 + 队列可合并性评估。

发现（解析前实测）：6 部芯片系列（GPL2112/6082/11648）series_matrix 内嵌完整表达表
（17k-45k 行）；2 部 RNA-seq（GPL15749: GSE131607/GSE75379）内嵌表为空壳（仅表头）
→ RNA-seq 轨需下载 counts/suppl 文件（下一步），本脚本先完成 6 部芯片轨解析。

输出：
  work/data/geo_parsed/<GSE>_expr.tsv.gz   （基因 x 样本 矩阵）
  work/output/de_queue_design.tsv          （样本级设计表：GSE/GSM/分组/组织/时间）
  work/output/de_queue_mergeability.json   （可合并性评估矩阵）
"""
import gzip, json, os, re, sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEO = os.path.join(HERE, 'data', 'geo')
PARSED = os.path.join(HERE, 'data', 'geo_parsed')
OUT = os.path.join(HERE, 'output')
os.makedirs(PARSED, exist_ok=True)

CHIPPED = ['GSE9685', 'GSE24560', 'GSE25319', 'GSE24217', 'GSE93082', 'GSE15025', 'GSE34031']

# 分组规则：按 char 字段值模式抽取 case/control（各系列字段名不同，逐一映射）
GROUP_RULES = {
    'GSE9685':   {'field': 'dose', 'case': lambda v: v != 'control', 'control': lambda v: v == 'control'},
    'GSE24560':  {'field': 'treatment', 'case': lambda v: v in ('escherichia coli', 'staphylococcus aureus'),
                  'control': lambda v: v == 'control'},
    'GSE25319':  {'field': 'infection', 'case': lambda v: v == 'post-infection', 'control': lambda v: v == 'pre-infection'},
    'GSE24217':  {'field': 'treatment group', 'case': lambda v: v.startswith('e.coli'), 'control': lambda v: v.startswith('control')},
    'GSE93082':  {'field': None, 'case': None, 'control': None},   # 田间：分组在 supplementary，标记 manual
    'GSE15025':  {'field': 'infection', 'case': lambda v: v.startswith('e.coli'), 'control': lambda v: v == 'uninfected'},
    'GSE34031':  {'field': 'disease state', 'case': lambda v: 'infected' in v, 'control': lambda v: 'healthy' in v},
}


def parse_matrix(gse):
    raw = gzip.open(os.path.join(GEO, f'{gse}_series_matrix.txt.gz'), 'rt',
                    encoding='utf-8', errors='replace').read()
    i = raw.find('series_matrix_table_begin')
    j = raw.find('series_matrix_table_end', i)
    lines = raw[i:j].strip().split('\n')
    # 首行可能是 'series_matrix_table_begin' 标记本身 → 找含 ID_REF 的一行作表头
    hdr_idx = 0
    for k, ln in enumerate(lines[:3]):
        if 'ID_REF' in ln:
            hdr_idx = k
            break
    header = [c.strip().strip('"') for c in lines[hdr_idx].rstrip('\r').split('\t')]
    samples = header[1:]
    ids, rows = [], []
    for ln in lines[hdr_idx + 1:]:
        parts = ln.rstrip('\r').split('\t')
        if not parts or parts[0].strip('"') in ('ID_REF', ''):
            continue
        ids.append(parts[0].strip('"'))
        try:
            vals = [np.nan if p.strip('"\r') in ('', 'null', 'NULL') else float(p)
                    for p in parts[1:len(samples) + 1]]
        except ValueError:
            vals = [np.nan] * len(samples)
        if len(vals) < len(samples):
            vals += [np.nan] * (len(samples) - len(vals))
        rows.append(vals[:len(samples)])
    # 样本级注解
    def col(field):
        m = re.search(rf'^!Sample_geo_accession\s+(.*)$', raw, re.M)
        gsm = [v.strip('"') for v in m.group(1).split('\t')] if m else []
        vals = {}
        for cm in re.finditer(rf'^!Sample_characteristics_ch1\s+(.*)$', raw, re.M):
            cs = [v.strip().strip('"') for v in cm.group(1).split('\t')]
            if any(field in c.lower() for c in cs if ':' in c):
                vals = {g: (c.split(':', 1)[1].strip() if ':' in c else c)
                        for g, c in zip(gsm, cs)}
                break
        return vals
    titles_m = re.search(r'^!Sample_title\s+(.*)$', raw, re.M)
    titles = [v.strip('"') for v in titles_m.group(1).split('	')] if titles_m else []
    titles = dict(zip([s for s in samples], titles))
    return {'samples': samples, 'ids': ids, 'mat': np.array(rows, dtype=np.float32),
            'titles': titles,
            'char': {f: col(f) for f in ('dose', 'treatment', 'infection', 'treatment group',
                                          'disease state', 'tissue', 'time point', 'time', 'time (h)')}}


def main():
    design_rows = []
    merge = {}
    for gse in CHIPPED:
        p = parse_matrix(gse)
        ids, samples, mat = p['ids'], p['samples'], p['mat']
        # 缺失值统计与过滤（保留 >=80% 样本非缺的探针）
        ok_frac = (~np.isnan(mat)).mean(axis=1)
        keep = ok_frac >= 0.6
        mat_f = mat[keep]
        # 缺失 = 探针中位数填充（GSE9685 空格率 43%，探针级中位稳健）
        med = np.nanmedian(mat_f, axis=1, keepdims=True)
        idx = np.where(np.isnan(mat_f))
        mat_f[idx] = np.take(med, idx[0])
        # z-score 平台内标准化（探针方向保留）
        mu = mat_f.mean(axis=1, keepdims=True)
        sd = mat_f.std(axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        z = (mat_f - mu) / sd
        out = os.path.join(PARSED, f'{gse}_expr.tsv.gz')
        with gzip.open(out, 'wt', encoding='utf-8') as f:
            f.write('probe\t' + '\t'.join(samples) + '\n')
            for pid, row in zip([ids[k] for k in range(len(ids)) if keep[k]], z):
                f.write(pid + '\t' + '\t'.join(f'{v:.4f}' for v in row) + '\n')
        # 分组
        rule = GROUP_RULES[gse]
        groups = {}
        if rule['field']:
            vals = p['char'].get(rule['field'], {})
            for s in samples:
                v = vals.get(s, '')
                v_l = str(v).lower()
                if rule['control'](v_l):
                    groups[s] = 'control'
                elif rule['case'](v_l):
                    groups[s] = 'case'
                else:
                    groups[s] = 'excluded'
        elif gse == 'GSE93082':
            for s in samples:
                tl = p['titles'].get(s, '').lower()
                if 'healthy' in tl:
                    groups[s] = 'control'
                elif 'mastitis' in tl:
                    groups[s] = 'case'
                else:
                    groups[s] = 'excluded'  # endometritis 非乳腺炎
        else:
            for s in samples:
                groups[s] = 'manual'
        n_case = sum(1 for v in groups.values() if v == 'case')
        n_ctl = sum(1 for v in groups.values() if v == 'control')
        for s in samples:
            design_rows.append({'gse': gse, 'gsm': s, 'group': groups[s],
                                'tissue_char': p['char'].get('tissue', {}).get(s, ''),
                                'time_char': (p['char'].get('time point', {}).get(s)
                                              or p['char'].get('time', {}).get(s)
                                              or p['char'].get('time (h)', {}).get(s, ''))})
        merge[gse] = {'probes_total': len(ids), 'probes_kept': int(keep.sum()),
                      'samples': len(samples), 'case': n_case, 'control': n_ctl,
                      'expr_gz': os.path.basename(out)}
        print(f"{gse}: probes {len(ids)}->{int(keep.sum())}, samples {len(samples)}, "
              f"case {n_case} / control {n_ctl}")
    # 设计表
    with open(os.path.join(OUT, 'de_queue_design.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv_writer = __import__('csv').DictWriter(f, fieldnames=['gse', 'gsm', 'group', 'tissue_char', 'time_char'], delimiter='\t')
        w.writeheader()
        w.writerows(design_rows)
    json.dump(merge, open(os.path.join(OUT, 'de_queue_mergeability.json'), 'w', encoding='utf-8'), indent=1)
    print('saved design + mergeability')


if __name__ == '__main__':
    main()
