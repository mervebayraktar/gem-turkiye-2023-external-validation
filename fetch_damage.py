"""Reproduce the observed Dataverse download without guessing filenames/schema."""
import hashlib
import json
import requests
from pipeline import RAW,ROOT,dump,fingerprint


def main():
    RAW.mkdir(parents=True,exist_ok=True)
    url='https://dataverse.harvard.edu/api/access/datafile/12080446'
    path=RAW/'rapor_final_v2.xlsx'
    if not path.exists():
        r=requests.get(url,timeout=120);r.raise_for_status();path.write_bytes(r.content)
    if hashlib.md5(path.read_bytes()).hexdigest()!='24957ee55a0b69ee9cd9cdf2a7325faa':
        raise ValueError('Published V1 file MD5 mismatch; do not analyse changed bytes silently')
    manifest=dict(doi='10.7910/DVN/OQCPWF',version='1.0',file_id=12080446,url=url,
                  sha256=fingerprint(path),md5='24957ee55a0b69ee9cd9cdf2a7325faa',
                  license='CC BY-NC-ND 4.0',license_evidence='https://api.datacite.org/dois/10.7910/DVN/OQCPWF')
    dump(ROOT/'outputs/tables/damage_source_manifest.json',manifest)
    # Cached full inventory was obtained with this exact public search query.
    inv=RAW/'dataverse_file_inventory.json'
    if not inv.exists():
        items=[];start=0
        while True:
            r=requests.get('https://dataverse.harvard.edu/api/search',params=dict(q='*',type='file',fq='parentId:12079754',
                              per_page=1000,start=start,sort='name',order='asc'),timeout=120)
            r.raise_for_status();data=r.json()['data'];items.extend(data['items'])
            start+=len(data['items'])
            if start>=data['total_count']:break
            if not data['items']:raise RuntimeError('Incomplete inventory pagination')
        if len({f['file_id'] for f in items})!=len(items):raise ValueError('Duplicate file ids in pagination')
        dump(inv,items)
    print('Dataverse V1 workbook verified; raw workbook unchanged')


if __name__=='__main__':main()
