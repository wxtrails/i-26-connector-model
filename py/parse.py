import xml.etree.ElementTree as ET, sys, json
ns={'k':'http://www.opengis.net/kml/2.2'}
t=ET.parse(sys.argv[1]); r=t.getroot()
out={}
for f in r.iter('{http://www.opengis.net/kml/2.2}Folder'):
    name=f.find('k:name',ns).text
    if name=='Levels': continue
    lines=[];polys=[]
    for ls in f.iter('{http://www.opengis.net/kml/2.2}LineString'):
        c=ls.find('k:coordinates',ns).text.split()
        lines.append([[float(v) for v in p.split(',')[:3]] for p in c])
    for pg in f.iter('{http://www.opengis.net/kml/2.2}Polygon'):
        rings=[]
        for lr in pg.iter('{http://www.opengis.net/kml/2.2}LinearRing'):
            c=lr.find('k:coordinates',ns).text.split()
            rings.append([[float(v) for v in p.split(',')[:3]] for p in c])
        polys.append(rings)
    out[name]={'lines':lines,'polys':polys}
    pts=[p for l in lines for p in l]+[p for pg in polys for rg in pg for p in rg]
    if pts:
        xs=[p[0] for p in pts]; ys=[p[1] for p in pts]; zs=[p[2] for p in pts]
        print(f"{name:45s} lines={len(lines):6d} polys={len(polys):6d} pts={len(pts):8d} lon[{min(xs):.4f},{max(xs):.4f}] lat[{min(ys):.4f},{max(ys):.4f}] z[{min(zs)},{max(zs)}]")
json.dump(out,open(sys.argv[2],'w'))
