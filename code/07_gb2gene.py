#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM2 Stage 2a-final: GPL2112 探针→基因三通道映射（符号/Entrez 直配 2,410 + GB_ACC 批量反查）。

GPL6082 同理（1,200 有 GB_LIST）。GB_ACC→基因走 NCBI efetch nucleotide→gene 或
efetch protein：批量 500/请求，绕代理。先建 GB_ACC 候选清单，efetch Batch Entrez。
简化方案：用 NCBI eutils efetch db=nucleotide & id=GB_ACC 列表 → 拿 /gene 或
/db_xref="GeneID:xxx"—— Affy 牛阵列探针对应 mRNA GB。
输出 work/data/platform/gb2gene.json (gb_acc -> {symbol, entrez})
"""
import gzip, json, os, re, sys, time, urllib.request, urllib.parse

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDIR = os.path.join(HERE, 'data', 'platform')
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
OPENER.addheaders = [('User-Agent', 'VM2-gate/1.0')]


def load_gbs(gpl, col='GB_ACC'):
    text = gzip.open(os.path.join(PDIR, 'soft', f'{gpl}_family.soft.gz'), 'rt', encoding='utf-8', errors='replace').read()
    cols, in_table = None, False
    gbs = {}
    for ln in text.split('\n'):
        if ln.startswith('!platform_table_begin'):
            in_table = True
        elif ln.startswith('!platform_table_end'):
            in_table = False
        elif in_table:
            if cols is None:
                cols = [c.strip() for c in ln.lstrip('#!').split('\t')]
            else:
                parts = ln.rstrip('\r').split('\t')
                row = dict(zip(cols, parts))
                pid = (row.get('ID') or row.get('D') or '').strip('"')
                gb = (row.get(col) or '').strip().strip('"').split('.')[0]
                if pid and gb:
                    gbs[pid] = gb
    return gbs


def efetch_geneids(gb_list):
    # efetch db=nucleotide rettype=gb -> 提取 /db_xref=GeneID 与 /gene= 标签
    out = {}
    B = 400
    for i in range(0, len(gb_list), B):
        chunk = gb_list[i:i + B]
        url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=text&id='
               + ','.join(chunk))
        for attempt in range(5):
            try:
                t = OPENER.open(url, timeout=120).read().decode('utf-8', errors='replace')
                break
            except Exception as e:
                print('  efetch retry', attempt + 1, type(e).__name__, flush=True)
                time.sleep(5 + attempt * 5)
        else:
            continue
        # 分记录（// 终止符；ACCESSION/VERSION 均可取登录号）
        for rec in t.split('\n//'):
            m_acc = re.search(r'(?:ACCESSION|VERSION)\s+(\S+)', rec)
            if not m_acc:
                continue
            acc = m_acc.group(1).split('.')[0]
            geneids = set(re.findall(r'/db_xref="GeneID:(\d+)"', rec))
            gene = re.search(r'/gene="([^"]+)"', rec)
            prod = re.search(r'/product="([^"]+)"', rec)
            if geneids:
                out[acc] = {'entrez': sorted(geneids)[0], 'gene': (gene.group(1) if gene else '')}
        print(f'  chunk {i // B + 1}/{(len(gb_list) - 1) // B + 1}: cum {len(out)} mapped', flush=True)
        time.sleep(1.2)
    return out


def main():
    cache = os.path.join(PDIR, 'gb2gene.json')
    done = json.load(open(cache)) if os.path.exists(cache) else {}
    for gpl in ('GPL2112', 'GPL6082'):
        gbs = load_gbs(gpl)
        uniq = sorted(set(gbs.values()))
        todo = [g for g in uniq if g not in done]
        print(f'{gpl}: {len(gbs)} probes, {len(uniq)} unique GB, {len(todo)} to fetch', flush=True)
        got = efetch_geneids(todo)
        done.update(got)
        json.dump(done, open(cache, 'w'))
        print(f'  cache now {len(done)} GB entries')


if __name__ == '__main__':
    main()
