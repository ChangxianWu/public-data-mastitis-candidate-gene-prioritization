#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 DE meta 队列重选与批量元数据采集（Stage 1）。

重选准则（mRNA 差异表达、乳腺相关组织/细胞、感染 vs 对照设计）：
  - 类型 = Expression profiling（池内 60 部全是，但含 miRNA/methylation 伪装 →
    用 series_type 细分 + 平台类型二次核验）
  - 组织 = mammary gland / mammary epithelial cells / milk cells / blood（宿主响应侧）
  - 设计 = 病原挑战/自然感染 vs 健康（DE meta 主线）；营养/能量平衡类排除
代表集（覆盖 3 组织 × 3 病原，n 样本 ≥ 9）：
  GSE9685(S aureus mammary tissue) GSE24560(E coli+S aureus pbMEC) GSE25319(E coli blood)
  GSE24217(E coli mammary) GSE93082(milk cells field) GSE131607(host response mammary)
  GSE15025(E coli mammary+blood) GSE34031(field mastitis mammary) GSE75379(E coli liver+udder)
输出: work/output/de_meta_queue_metadata.json + work/data/geo/<GSE>_series_matrix.txt.gz
"""
import gzip, json, os, re, sys, time
os.environ["PYTHONUNBUFFERED"] = "1"
import urllib.request

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, 'output')
DATA = os.path.join(HERE, 'data', 'geo')
os.makedirs(DATA, exist_ok=True)

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
OPENER.addheaders = [('User-Agent', 'VM2-gate/1.0 (DE-meta queue selection)')]

QUEUE = ['GSE9685', 'GSE24560', 'GSE25319', 'GSE24217', 'GSE93082',
         'GSE131607', 'GSE15025', 'GSE34031', 'GSE75379']


def fetch(url, tries=5, timeout=120):
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
            print('    retry', i + 1, type(e).__name__)
            time.sleep(4 + i * 5)
    raise last


def gse_dir(gse):
    # GEO 官方规则: 目录 = GSE + str(n // 1000) + 'nnn'（GSE9685->GSE9nnn, 24560->GSE24nnn, 131607->GSE131nnn）
    n = int(re.match(r'GSE(\d+)', gse).group(1))
    return 'GSE{}nnn'.format(n // 1000)


def parse_matrix(gse):
    dst = os.path.join(DATA, f'{gse}_series_matrix.txt.gz')
    if not os.path.exists(dst) or os.path.getsize(dst) < 500:
        url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_dir(gse)}/{gse}/matrix/{gse}_series_matrix.txt.gz'
        b = fetch(url)
        open(dst, 'wb').write(b)
        time.sleep(1.2)
    raw = gzip.open(dst, 'rt', encoding='utf-8', errors='replace').read()
    rec = {'gse': gse, 'bytes': os.path.getsize(dst)}
    for key in ('!Series_title', '!Series_type', '!Series_platform_id', '!Series_sample_taxid'):
        m = re.search(rf'^{re.escape(key)}\t"?(.*?)"?\s*$', raw, re.M)
        if m:
            rec[key] = m.group(1).split('\t')[0].strip('"')[:180]
    ids = re.search(r'!Sample_geo_accession\s+(.*)', raw)
    rec['n_samples'] = len(ids.group(1).split('\t')) if ids else 0
    titles = re.search(r'!Sample_title\s+(.*)', raw)
    rec['sample_title_head'] = [v.strip('"') for v in titles.group(1).split('\t')[1:6]] if titles else []
    chars = re.findall(r'!Sample_characteristics_ch1\s+(.*)', raw)
    fields = []
    for line in chars:
        vals = [v.strip().strip('"') for v in line.split('\t')[1:] if v.strip()]
        if vals:
            uniq = sorted({v.lower() for v in vals})
            fields.append(uniq[:8])
    rec['char_fields'] = fields
    rec['has_expr_table'] = 'series_matrix_table_begin' in raw
    rec['suppl_expression_hint'] = bool(re.search(r'!Series_supplementary_file.*?_raw|FPKM|counts|tpm', raw, re.I))
    return rec


def main():
    out, failed = {}, []
    for gse in QUEUE:
        try:
            rec = parse_matrix(gse)
        except Exception as e:
            failed.append({'gse': gse, 'error': f'{type(e).__name__}: {str(e)[:80]}'})
            print(f'{gse}: FAILED {type(e).__name__}, continuing')
            time.sleep(2)
            continue
        out[gse] = rec
        print(f"{gse}: n={rec['n_samples']} platform={rec.get('!Series_platform_id')} "
              f"type={rec.get('!Series_type', '')[:42]} expr_table={rec['has_expr_table']}", flush=True)
        for f in rec['char_fields']:
            print('    char:', f[:6])
        json.dump({'queue': out, 'failed': failed},
                  open(os.path.join(OUT, 'de_meta_queue_metadata.json'), 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=1)
        time.sleep(1.0)
    json.dump({'queue': out, 'failed': failed},
              open(os.path.join(OUT, 'de_meta_queue_metadata.json'), 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1)
    print(f'saved: {len(out)} ok / {len(failed)} failed')


if __name__ == '__main__':
    main()
