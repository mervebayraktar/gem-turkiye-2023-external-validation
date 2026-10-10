"""Download pinned public GEM/USGS sources; preserve URLs and SHA256.
No damage-data URL/file schema is invented. Network failures remain explicit.
"""
import argparse
import json
from datetime import datetime, timezone
import requests
from pipeline import RAW, ROOT, SHA, fingerprint, dump

EVENTS = {'event1': ('us6000jllz', '1756921940993'),
          'event2': ('us6000jlqa', '1756575631263'),
          'event3': ('us6000jlm1', '1679026302194'),
          'event4': ('us6000jqcn', '1728481576510')}


def acquire():
    RAW.mkdir(parents=True,exist_ok=True)
    urls={}
    for c in ['structural','nonstructural','contents','fatalities']:
        urls[f'gem_vulnerability_{c}.xml']=f'https://raw.githubusercontent.com/gem/global_vulnerability_model/{SHA}/Europe/Turkey/vulnerability_{c}.xml'
    for f in ['README.md','LICENSE.txt']:
        urls['gem_'+f]=f'https://raw.githubusercontent.com/gem/global_vulnerability_model/{SHA}/{f}'
    urls['tree.json']=f'https://api.github.com/repos/gem/global_vulnerability_model/git/trees/{SHA}?recursive=1'
    for label,(event,update) in EVENTS.items():
        # Snapshot is pinned. Event detail endpoint itself can change, and is labelled as such.
        for name in ['grid.xml','uncertainty.xml','info.json','rupture.json']:
            urls[f'{label}_{name}']=f'https://earthquake.usgs.gov/product/shakemap/{event}/us/{update}/download/{name}'
        urls[f'{label}.json']=f'https://earthquake.usgs.gov/earthquakes/feed/v1.0/detail/{event}.geojson'
    urls['datacite.txt']='https://api.datacite.org/dois/10.7910/DVN/OQCPWF'
    manifest=[]
    for name,url in urls.items():
        path=RAW/name
        record=dict(file=str(path.relative_to(ROOT)),url=url,checked_utc=datetime.now(timezone.utc).isoformat())
        try:
            if not path.exists():
                r=requests.get(url,timeout=120);r.raise_for_status();path.write_bytes(r.content)
                record['action']='downloaded'
            else:record['action']='existing_bytes_preserved'
            record.update(sha256=fingerprint(path),bytes=path.stat().st_size,status='available')
        except requests.RequestException as e:
            record.update(status='bulunamadı',error=str(e))
        manifest.append(record)
    dump(ROOT/'outputs/tables/source_manifest.json',manifest)
    if any(r['status']!='available' for r in manifest):
        raise SystemExit('Bazı kaynaklar bulunamadı; manifesti inceleyin')
    print(f'{len(manifest)} sources recorded; existing bytes preserved')


if __name__=='__main__':acquire()
