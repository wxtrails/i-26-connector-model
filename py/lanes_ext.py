"""Carry the plan's lane layout onto the roads that continue past the plan limits.

Where the plan stops, the model continues on simplified roads ("extensions", see ext.py), and there is often a stretch
of pavement between the last drawn lane line and the extension. Here we read a lane template (lane line, edge and
centre positions across the road) from the nearest clean cross-section of the plan, then use it to
  * size the extension pavement (instead of a fixed width),
  * paint lane lines along the gap stretch and along the extension, tapering to the extension's default width.
"""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg, shapely.ops as so
from ext import EXT
SKIP={'south'}                      # the I-40 end uses the 2018 KMZ (no lane lines); handled elsewhere
DEFAULT_HW={'i240e':7.6,'north':7.6,'pate':9.2,'patw':9.2}
TAPER=(220.0,300.0)                 # template held for the first 220 m of the extension, then tapers to the default width
def _lines(L,key):
    out=[]
    for it in L.get(key,[]):
        a=np.asarray(it[1] if isinstance(it,tuple) else it,float)
        if len(a)>=2 and np.hypot(*np.diff(a,axis=0).T).sum()>0.5: out.append(sg.LineString(a))
    return out
def _frame(P):
    t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; return t,np.c_[-t[:,1],t[:,0]]   # normal = left of travel
def _resample(P,ds=2.0):
    ls=sg.LineString(P); s=np.arange(0,ls.length,ds); s=np.r_[s,ls.length]
    return s,np.array([ls.interpolate(x).coords[0] for x in s])
def plan_side(chain_xy,joint):
    """plan alignment as a polyline that starts at the joint and runs into the plan"""
    P=np.asarray(chain_xy,float)
    if np.hypot(*(P[-1]-joint))<np.hypot(*(P[0]-joint)): P=P[::-1]
    return P
class Pts:
    """lines as dense points, so a cross-section also catches dashed lines whose dashes miss it"""
    def __init__(self,lines,step=1.0):
        P=[]
        for g in lines:
            n=max(2,int(g.length/step)+1); P.append(np.array([g.interpolate(x).coords[0] for x in np.linspace(0,g.length,n)]))
        self.P=np.vstack(P) if P else np.zeros((0,2))
        from scipy.spatial import cKDTree; self.tree=cKDTree(self.P)
def cross(pts,p,t,n,half,win=7.0):
    idx=pts.tree.query_ball_point(p,np.hypot(half,win))
    if not idx: return []
    Q=pts.P[idx]-p; al=Q@t; la=Q@n; m=(np.abs(al)<=win)&(np.abs(la)<=half)
    la=np.sort(la[m]); out=[]
    for o in la:                      # cluster offsets within 0.6 m
        if out and o-out[-1][-1]<0.6: out[-1].append(o)
        else: out.append([o])
    return [float(np.mean(c)) for c in out]
def find_template(Pplan,LL,EO,Graw,hw_plan,maxd=260.0):
    """walk into the plan until a cross-section shows a stable lane layout; returns offsets relative to the plan alignment"""
    lt=Pts(LL); et=Pts(EO)
    s,P=_resample(Pplan,4.0); T,N=_frame(P); half=hw_plan+10
    prev=None
    for i in range(len(s)):
        if s[i]>maxd: break
        ln=cross(lt,P[i],T[i],N[i],half); eo=cross(et,P[i],T[i],N[i],half,win=2.0)
        if len(ln)<2: prev=None; continue
        lo=[e for e in eo if e<min(ln)-1.0]; hi=[e for e in eo if e>max(ln)+1.0]
        if not lo or not hi: prev=None; continue
        tpl={'lines':ln,'edges':(max(lo),min(hi)),'s':float(s[i]),'p':P[i]}
        if prev and len(prev['lines'])==len(ln) and abs(prev['edges'][0]-tpl['edges'][0])<1.0 and abs(prev['edges'][1]-tpl['edges'][1])<1.0:
            return prev
        prev=tpl
    return None
def gap_end(Pplan,LL,Graw,tpl,maxd=200.0):
    """distance from the joint into the plan over which no lane line is drawn and the road is not an intersection"""
    lt=Pts(LL); s,P=_resample(Pplan,2.0); T,N=_frame(P); lo,hi=tpl['edges']; w=hi-lo
    for i in range(len(s)):
        if s[i]>maxd: return s[i]
        ln=[o for o in cross(lt,P[i],T[i],N[i],max(abs(lo),abs(hi))+1,win=1.0) if lo-0.5<o<hi+0.5]
        if ln: return s[i]
        seg=sg.LineString([P[i]+N[i]*(lo-12),P[i]+N[i]*(hi+12)]); dr=seg.intersection(Graw).length
        if dr>w+8: return s[i]           # pavement much wider than the road: an intersection, stop before it
    return s[-1]
def classify(tpl):
    """kinds per template line: 'y' solid yellow beside the centre, 'w' dashed white elsewhere"""
    ln=np.array(tpl['lines']); kinds=['w']*len(ln)
    if len(ln)>=2:
        g=np.diff(ln); k=int(np.argmax(g))
        if g[k]>4.5:                    # median, island or turn lane: bound it with yellow
            kinds[k]='y'; kinds[k+1]='y'
        else:
            j=int(np.argmin(np.abs(ln-np.mean(tpl['edges'])))); kinds[j]='yy'
    return list(zip(ln.tolist(),kinds))
def build(chains=None,E=None,verbose=True,hwf=None):
    if chains is None: chains=pickle.load(open(W+'/profiles.pkl','rb'))['chains']
    if E is None: E=pickle.load(open(W+'/ext.pkl','rb'))
    L=pickle.load(open(W+'/pdf_layers.pkl','rb'))
    Graw=pickle.load(open(W+'/sheets_plan.pkl','rb'))['G']
    LL=_lines(L,('Prop Lane Lines','s')); EO=_lines(L,('Prop EOT','s'))
    if hwf is None: from profiles2 import hw as hwf
    HW=hwf
    out={}
    for k,v in E.items():
        if k in SKIP: continue
        joint=np.array(EXT[k][0],float); best=None
        for c,d in chains.items():
            if not c.startswith('b'): continue
            for e in (0,-1):
                dd=np.hypot(*(np.asarray(d['xy'])[e]-joint))
                if dd<20 and (best is None or dd<best[0]): best=(dd,c)
        if best is None: print(k,'no plan alignment at the joint'); continue
        c=best[1]; Pp=plan_side(chains[c]['xy'],joint)
        tpl=find_template(Pp,LL,EO,Graw,HW(c,chains))
        if tpl is None: print(k,'no lane template found in',c); continue
        ge=gap_end(Pp,LL,Graw,tpl)
        refs=[v['out']] if not v['dual'] else [v['out'],v['in'][::-1]]
        # lateral position of each carriageway reference at the joint, relative to the plan alignment (normal = left looking outward)
        t0=joint-Pp[min(3,len(Pp)-1)]; t0/=np.hypot(*t0); n0=np.array([-t0[1],t0[0]])
        # template normal points left looking INTO the plan; flip to "left looking outward"
        flip=-1.0
        lines=[(flip*o,kd) for o,kd in classify(tpl)]; e=sorted([flip*tpl['edges'][0],flip*tpl['edges'][1]])
        pos=[float((np.asarray(r[0])-joint)@n0) for r in refs]
        out[k]={'plan':c,'tpl_s':tpl['s'],'gap':ge,'lines':sorted(lines),'edges':tuple(e),'refs':[np.asarray(r,float) for r in refs],'pos':pos,'dual':v['dual'],
                'hw_def':DEFAULT_HW[k],'Pplan':Pp}
        if verbose: print(k,'plan',c,'template at %.0f m'%tpl['s'],'gap %.0f m'%ge,'edges',np.round(e,1),'lines',[(round(o,1),kd) for o,kd in sorted(lines)],'refs at',np.round(pos,1))
    return out
def _smooth(x):
    x=np.clip(x,0,1); return x*x*(3-2*x)
def groups(t,ds=2.0):
    """one painted group per carriageway: outward path C (gap stretch + extension), and per point the lateral offsets"""
    Pp=t['Pplan']; sp,Pq=_resample(Pp,ds); gap=Pq[sp<=t['gap']][::-1]          # from the gap end back to the joint
    lines=t['lines']; lo,hi=t['edges']; out=[]
    def join(gp,R):
        # drop gap points that reach past the start of the extension, so the path never doubles back
        if len(gp)==0 or len(R)<2: return gp
        tR=R[min(3,len(R)-1)]-R[0]; tR/=np.hypot(*tR)+1e-9
        return gp[(gp-R[0])@tR<-2.0]
    median=sum(1 for o,kd in lines if kd=='y')>=2
    if t['dual'] and not median:          # two centrelines but the plan road has no median: treat as one road on their midline
        a,b=t['refs']; f=np.linspace(0,1,200)
        la,lb=sg.LineString(a),sg.LineString(b)
        mid=np.array([(np.array(la.interpolate(x,normalized=True).coords[0])+np.array(lb.interpolate(x,normalized=True).coords[0]))/2 for x in f])
        t=dict(t,dual=False,refs=[mid],pos=[float(np.mean(t['pos']))])
    if not t['dual']:
        ref=t['refs'][0]; sr,R=_resample(ref,ds); _,Nr=_frame(R); w=1-_smooth(sr/150.0)
        R=R-Nr*(t['pos'][0]*w)[:,None]                                        # bend the extension onto the plan alignment at the joint
        mid=(lo+hi)/2; gp=join(gap,R); C=np.vstack([gp,R]); s=np.r_[np.zeros(len(gp)),sr]
        rec=np.r_[np.zeros(len(gp)),mid*(1-w)]                               # recentre the template on the extension's own centreline
        grp=[(C,s,[(o,kd) for o,kd in lines],lo,hi,rec)]
    else:
        ym=[o for o,kd in lines if kd.startswith('y')]; cmed=np.mean(ym) if ym else np.mean(t['pos'])
        grp=[]
        for j,ref in enumerate(t['refs']):
            pj=t['pos'][j]; side=np.sign(pj-cmed) or 1.0
            mine=[(o-pj,kd) for o,kd in lines if np.sign(o-cmed)==side]
            ys=[o-pj for o,kd in lines if kd.startswith('y') and np.sign(o-cmed)==side]
            outer=(hi if side>0 else lo)-pj; inner=(ys[0]-side*1.0) if ys else (cmed-pj)
            a,b=sorted([outer,inner])
            sr,R=_resample(ref,ds); gp=join(gap+_frame(gap)[1]*pj if len(gap)>1 else gap,R); C=np.vstack([gp,R]); s=np.r_[np.zeros(len(gp)),sr]
            grp.append((C,s,mine,a,b,np.zeros(len(s))))
    for C,s,ln,a,b,rec in grp:
        hw=t['hw_def']; f=_smooth((s-TAPER[0])/(TAPER[1]-TAPER[0]))
        c0=(a+b)/2
        A=(a-rec)*(1-f)+(-hw)*f if not t['dual'] else a*(1-f)+(c0-hw)*f
        B=(b-rec)*(1-f)+(hw)*f if not t['dual'] else b*(1-f)+(c0+hw)*f
        out.append({'C':C,'s':s,'A':A,'B':B,'lines':[(np.full(len(s),o)-rec,kd) for o,kd in ln]})
    return out
def paint(T):
    """polylines to paint: list of (Nx2, kind) with kind 'dash' (white), 'white', 'yellow'"""
    res=[]
    def runs(mask):
        idx=np.where(mask)[0]
        return [r for r in np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(r)>=3] if len(idx) else []
    for k,t in T.items():
        for g in groups(t):
            C,s,A,B=g['C'],g['s'],g['A'],g['B']; _,N=_frame(C)
            for side in (A+0.3,B-0.3):
                res.append((C+N*side[:,None],'white'))
            for o,kd in g['lines']:
                ok=(o>A+0.8)&(o<B-0.8)
                for r in runs(ok):
                    P=C[r]+N[r]*o[r][:,None]
                    if kd=='w': res.append((P,'dash'))
                    elif kd=='y': res.append((P,'yellow'))
                    else:
                        res.append((C[r]+N[r]*(o[r]-0.15)[:,None],'yellow')); res.append((C[r]+N[r]*(o[r]+0.15)[:,None],'yellow'))
    return res
def pavement(T):
    """pavement polygons for the gap stretch and the extension, from the template edges"""
    polys=[]
    for k,t in T.items():
        for g in groups(t):
            C,A,B=g['C'],g['A'],g['B']; _,N=_frame(C)
            ring=np.vstack([C+N*(A-0.5)[:,None],(C+N*(B+0.5)[:,None])[::-1]])
            p=sg.Polygon(ring).buffer(0)
            if not p.is_empty: polys.append(p)
    return shapely.unary_union(polys) if polys else sg.Polygon()
if __name__=='__main__':
    T=build()
    pickle.dump(T,open(W+'/lanes_ext.pkl','wb'))
def stub_extensions(Graw,Gfinal,lines,maxlen=120.0,step=1.0,exclude=None,minlen=15.0):
    """Plan lane lines that stop where the drawn pavement stops, while the model's pavement carries on (filler pavement
    added past the plan limit): continue them straight along their end direction until the pavement ends."""
    shapely.prepare(Graw); shapely.prepare(Gfinal); Gb=Graw.boundary; out=[]
    for g in lines:
        a=np.asarray(g.coords)
        if len(a)<2 or g.length<1.0: continue
        for end in (0,-1):
            p=a[end]; q=a[1] if end==0 else a[-2]; t=p-q; L=np.hypot(*t)
            if L<1e-6: continue
            t=t/L
            if Gb.distance(sg.Point(*p))>3.0: continue
            ahead=p+t*3.0
            if Graw.contains(sg.Point(*ahead)) or not Gfinal.contains(sg.Point(*ahead)): continue
            n=0; pts=[p]
            while n*step<maxlen:
                n+=1; x=p+t*n*step
                if not Gfinal.contains(sg.Point(*x)) or Graw.contains(sg.Point(*x)): break
                if exclude is not None and exclude.contains(sg.Point(*x)): break
                pts.append(x)
            pts=np.array(pts[:-1]) if len(pts)>1 else np.array(pts)
            # only real continuations past the plan limit: long enough and leading well away from the drawn pavement
            if len(pts)>=2 and np.hypot(*(pts[-1]-pts[0]))>=minlen and Graw.distance(sg.Point(*pts[-1]))>=8.0: out.append(pts)
    return out
