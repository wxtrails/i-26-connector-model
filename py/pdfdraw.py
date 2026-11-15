import fitz, sys, pickle, collections
src,out=sys.argv[1],sys.argv[2]
doc=fitz.open(src); pg=doc[0]
print(pg.rect, pg.rotation)
dr=pg.get_drawings()
print(len(dr))
simp=[]
cnt=collections.Counter()
for d in dr:
    items=[]
    for it in d['items']:
        if it[0]=='l': items.append(('l',(it[1].x,it[1].y),(it[2].x,it[2].y)))
        elif it[0]=='c': items.append(('c',(it[1].x,it[1].y),(it[2].x,it[2].y),(it[3].x,it[3].y),(it[4].x,it[4].y)))
        elif it[0]=='re': r=it[1]; items.append(('re',(r.x0,r.y0),(r.x1,r.y1)))
        elif it[0]=='qu': q=it[1]; items.append(('qu',(q.ul.x,q.ul.y),(q.ur.x,q.ur.y),(q.lr.x,q.lr.y),(q.ll.x,q.ll.y)))
    key=(d['type'],tuple(round(c,2) for c in d['fill']) if d.get('fill') else None,tuple(round(c,2) for c in d['color']) if d.get('color') else None,round(d.get('width') or 0,2), d.get('dashes'))
    cnt[key]+=len(items)
    simp.append({'type':d['type'],'fill':d.get('fill'),'color':d.get('color'),'width':d.get('width'),'dashes':d.get('dashes'),'items':items,'rect':tuple(d['rect']),'fo':d.get('fill_opacity'),'layer':d.get('layer')})
pickle.dump(simp,open(out,'wb'))
for k,v in cnt.most_common(60): print(v,k)
print(set(s['layer'] for s in simp))
