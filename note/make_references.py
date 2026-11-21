"""Builds note/references.toml from the PubMed esummary records saved in note/sources/*.json (fetched from NCBI
E-utilities; each fetch is in the workbench record). No field is typed by hand: authors, title, journal, year,
volume, issue, pages and DOI come from the record. KEYS maps the note's citation key to its PMID. Software without a
PubMed record is listed in SOFTWARE. Run: python3 note/make_references.py"""
import json, glob, re
KEYS = {  # citation key -> PMID (order is irrelevant: build_pdf numbers by first citation)
 'lefebvre1995': '7813012', 'roy1995': '7813013', 'burglen1997': '8981949',
 'brzustowicz1990': '2320125', 'melki1990': '1970420', 'gilliam1990': '1972783',
 'lorson1999': '10339583', 'monani1999': '10369862', 'mercuri2022': '35927425',
 'finkel2017': '29091570', 'mendell2017': '29091557', 'prior2008': '18941424',
 'hendrickson2009': '19625283', 'sugarman2012': '21811307', 'feng2017': '28125085', 'luo2014': '23788250',
 'chen2020': '32066871', 'chen2023': '36669496', 'milligan2021': '33798739', 'mikhalchuk2025': '40076610',
 'campbell1997': '9245983', 'vandersteege1996': '8808598', 'wadman2020': '32954327', 'vorster2020': '32117462',
 'zwartkruis2025': '40119448', 'berzal2026': '42162012', 'suzuki2026': '42203797', 'logsdon2025': '40702183',
 'liao2023': '37165242', 'hickey2024': '37165083', 'chin2023': '37365340', 'yoo2025': '40205052',
 'rochette2001': '11354640', 'katoh2013': '23329690', 'minh2020': '32011700', 'hoang2018': '29077904',
 'virtanen2020': '32015543', 'mantel1967': '6018555',
 'bolognini2024': '39232174', 'yilmaz2024': '39418342', 'rozen2003': '12815433', 'bonito2021': '34244762',
 'harris2020': '32939066',
}
SOFTWARE = [dict(key='matplotlib', authors=['The Matplotlib Development Team'], title='Matplotlib: visualization with Python',
                 year=2025, url='https://matplotlib.org', cited='2026-11-14')]
HAND = [  # no PubMed or Europe PMC record: fields exactly as the researcher supplied them (2026-11-14), not from a record
 dict(key='anderson2001', authors=['Anderson MJ'], title='A new method for non-parametric multivariate analysis of variance',
      journal='Austral Ecology', year=2001, volume='26', pages='32--46', doi='10.1111/j.1442-9993.2001.01070.pp.x'),
 # sokal1958 not in the numbered list: build_pdf requires a doi or pmid for an article and the 1958 paper has neither.
 # Researcher's choice (2026-11-14): cite it as a plain in-text reference, 'Sokal and Michener, 1958', in Fig. 2's caption and Methods.
 # (Sokal RR, Michener CD. A statistical method for evaluating systematic relationships. Univ Kansas Sci Bull 38:1409-1438, as the researcher supplied it.)
 # dict(key='sokal1958', authors=['Sokal', 'Michener'], title='A statistical method for evaluating systematic relationships',
#     journal='University of Kansas Science Bulletin', year=1958, volume='38', pages='1409--1438'),
]
PREPRINTS = {  # citation key -> saved Europe PMC search record (preprints have no PubMed record)
 'lucas2026': 'note/sources/europepmc_hprc2_lucas2026.json',
}
rec = {}
for f in sorted(glob.glob('note/sources/*.json')):
    try: d = json.load(open(f))
    except Exception: continue
    r = d.get('result') if isinstance(d, dict) else None
    if not r or 'uids' not in r: continue
    for u in r['uids']: rec.setdefault(u, (r[u], f))
def pages(p):
    m = re.match(r'^(\d+)-(\d+)$', p or '')
    if not m: return p or None
    a, b = m.groups()
    if len(b) < len(a): b = a[:len(a) - len(b)] + b
    return f'{a}--{b}'
def q(s): return '"' + str(s).replace('\\', '\\\\').replace('"', '\\"') + '"'
out = ['# written by note/make_references.py from note/sources/*.json (PubMed esummary); do not edit']
missing = [k for k, p in KEYS.items() if p not in rec]
assert not missing, f'no saved PubMed record for {missing}'
for k, p in KEYS.items():
    r, f = rec[p]
    doi = next((a['value'] for a in r.get('articleids', []) if a['idtype'] == 'doi'), None)
    au = [a['name'] for a in r.get('authors', []) if a.get('authtype') == 'Author']
    year = int(re.match(r'\d{4}', r['pubdate']).group(0))
    t = r['title'].rstrip('.')
    out += ['', '[[ref]]', f'key = {q(k)}', 'kind = "article"', 'authors = [' + ', '.join(q(a) for a in au) + ']',
            f'title = {q(t)}', f'journal = {q(r["source"])}', f'year = {year}']
    for fld, v in [('volume', r.get('volume')), ('issue', (r.get('issue') or '').split('_')[0]),  # '2_Part_1' -> '2' (TeX-safe)
                    ('pages', pages(r.get('pages'))), ('doi', doi)]:
        if v: out.append(f'{fld} = {q(v)}')
    out += [f'pmid = {q(p)}', f'checked = {q("PubMed esummary, " + f)}']
for k, f in PREPRINTS.items():
    r = json.load(open(f))['resultList']['result'][0]
    assert 'Preprint' in r['pubTypeList']['pubType'], f
    au = [a['fullName'] for a in r['authorList']['author'] if 'fullName' in a]
    out += ['', '[[ref]]', f'key = {q(k)}', 'kind = "article"', 'authors = [' + ', '.join(q(a) for a in au) + ']',
            f'title = {q(r["title"].rstrip("."))}', f'journal = {q(r["bookOrReportDetails"]["publisher"] + " (preprint)")}',
            f'year = {int(r["pubYear"])}', f'doi = {q(r["doi"])}', f'checked = {q("Europe PMC record " + r["id"] + ", " + f)}']
for h in HAND:
    out += ['', '[[ref]]', f'key = {q(h["key"])}', 'kind = "article"', 'authors = [' + ', '.join(q(a) for a in h['authors']) + ']',
            f'title = {q(h["title"])}', f'journal = {q(h["journal"])}', f'year = {h["year"]}']
    out += [f'{k} = {q(h[k])}' for k in ('volume', 'pages', 'doi') if h.get(k)]
    out += [f'checked = {q("entered by hand from the researcher, 2026-11-14; no PubMed record")}']
for s in SOFTWARE:
    out += ['', '[[ref]]', f'key = {q(s["key"])}', 'kind = "software"', 'authors = [' + ', '.join(q(a) for a in s['authors']) + ']',
            f'title = {q(s["title"])}', f'year = {s["year"]}', f'url = {q(s["url"])}', f'cited = {q(s["cited"])}']
open('note/references.toml', 'w').write('\n'.join(out) + '\n')
print(len(KEYS) + len(PREPRINTS) + len(HAND) + len(SOFTWARE), 'references')
