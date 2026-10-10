"""Distance to the shallowest rupture-edge polyline, not Rrup or exact Rjb.

Uses the observed USGS MultiPolygon JSON. Local azimuthal-equidistant projection
centred on the event; planar polyline approximation in km is explicitly labelled.
"""
import json
import numpy as np
from pyproj import CRS, Transformer
from shapely.geometry import LineString, Point
from shapely.ops import unary_union


def trace_distance(path, lat, lon):
    data=json.loads(path.read_text(encoding='utf-8'))
    meta=data['metadata']
    crs=CRS.from_proj4(f"+proj=aeqd +lat_0={meta['lat']} +lon_0={meta['lon']} +datum=WGS84 +units=m")
    project=Transformer.from_crs('EPSG:4326',crs,always_xy=True)
    segments=[]
    for f in data['features']:
        geom=f['geometry']
        if geom['type']!='MultiPolygon':raise ValueError('Rupture geometry changed; review needed')
        for polygon in geom['coordinates']:
            for ring in polygon:
                a=np.asarray(ring,float)
                top=a[np.isclose(a[:,2],a[:,2].min()),:2]
                # Closing vertex would create an artificial diagonal shortcut.
                if len(top)>2 and np.allclose(top[-1],top[0]):top=top[:-1]
                if len(top)<2:continue
                x,y=project.transform(top[:,0],top[:,1])
                segments.append(LineString(np.column_stack([x,y])))
    if not segments:raise ValueError('Shallow rupture edge bulunamadı')
    trace=unary_union(segments)
    x,y=project.transform(np.asarray(lon).tolist(),np.asarray(lat).tolist())
    return np.array([Point(xx,yy).distance(trace)/1000 if np.isfinite(xx+yy) else np.nan for xx,yy in zip(x,y)])
