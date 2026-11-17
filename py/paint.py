"""US pavement marking roles for the plan's lane lines and edge-of-travel lines.

Standard US markings (MUTCD), as in the NCDOT typical sections: a solid yellow line on the left edge of each one-way
roadway (the median side of a divided road, the left side of a ramp), a solid white line on the right edge, dashed
white lines between lanes that run the same way, and a double yellow centre line between lanes that run opposite ways.

The plans draw the lines but not their role. Which way the lanes run comes from the roads' alignments and the plan's
lane arrows: an alignment whose arrows all point one way is a one-way roadway (a ramp, or one direction of a road
drawn as two alignments); otherwise it carries two-way traffic, the lanes on its right running its way. An edge line's
lanes lie on the side where the next line is a lane's width away."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
from scipy.spatial import cKDTree
LANE=(2.7,5.8)        # a lane's width (m): the distance from an edge line to the next line on its lanes' side
WIN=250.0             # arrows within this distance along an alignment set its traffic there

def _frame(P):
    t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9
    return t,np.c_[-t[:,1],t[:,0]]

class Traffic:
    """which way traffic runs at a point, from the alignment it belongs to"""
    def __init__(self,cls,chains,hw,arrows):
        P=[];T=[];C=[];H=[];I=[];self.mode={};self.both={}
        A=arrows
        for c,v in cls.items():
            Q=v['P']; t,nv=_frame(Q); h=hw(c,chains); s=v['s']
            # local traffic per station from the arrows on this alignment's carriageways
            fw=np.zeros(len(Q)); bw=np.zeros(len(Q)); rf=np.zeros(len(Q)); lb=np.zeros(len(Q)); nf=nb_=0
            if len(A):
                d,j=cKDTree(Q).query(A[:,:2])
                lat=np.einsum('ij,ij->i',A[:,:2]-Q[j],nv[j]); dot=np.einsum('ij,ij->i',A[:,2:4],t[j])
                m=(np.abs(lat)<=h+1)&(np.abs(dot)>0.85)&(d<h+6)
                nf=int((m&(dot>0)).sum()); nb_=int((m&(dot<0)).sum())
                for k,jj in zip(np.where(m)[0],j[m]):
                    w=np.abs(s-s[jj])<=WIN
                    if dot[k]>0: fw[w]+=1
                    else: bw[w]+=1
                    if dot[k]>0 and lat[k]<0: rf[w]+=1
                    if dot[k]<0 and lat[k]>0: lb[w]+=1
            tot=fw+bw
            md=np.zeros(len(Q),int)                    # 0 two-way, +1 one-way along, -1 one-way against
            one=tot>0
            # two-way needs arrows running its way on its right AND against it on its left; otherwise the majority
            # (a few arrows of a neighbouring roadway within reach do not make a one-way road two-way)
            two=(rf>0)&(lb>0)&(np.minimum(rf,lb)>=0.15*tot)
            md[one&~two&(fw>bw)]=1; md[one&~two&(bw>fw)]=-1
            # stations with no arrows nearby take the alignment's overall verdict
            if (~one).any():
                F,B=nf,nb_
                ov=1 if (F>0 and B<=0.15*F) else -1 if (B>0 and F<=0.15*B) else 0
                md[~one]=ov
            self.mode[c]=md; self.both[c]=(fw>0)&(bw>0)
            P.append(Q);T.append(t);C+= [c]*len(Q);H.append(np.full(len(Q),h));I.append(np.arange(len(Q)))
        self.P=np.vstack(P);self.T=np.vstack(T);self.C=np.array(C);self.H=np.concatenate(H);self.I=np.concatenate(I)
        self.tree=cKDTree(self.P)
    def owner(self,p,tl):
        """(chain, station, lateral offset +left, mode) of the alignment whose carriageway p lies in, parallel to tl"""
        idx=self.tree.query_ball_point(p,40.0)
        if not idx: return None
        idx=np.array(idx); dv=p-self.P[idx]; c=np.abs(self.T[idx]@tl)
        lat=dv[:,0]*(-self.T[idx,1])+dv[:,1]*self.T[idx,0]; along=np.abs(np.einsum('ij,ij->i',dv,self.T[idx]))
        ok=(c>0.85)&(np.abs(lat)<=self.H[idx]+3.0)&(along<4.0)
        if not ok.any(): return None
        idx=idx[ok]; lat=lat[ok]; c=c[ok]
        k=np.argmin(np.maximum(np.abs(lat)-self.H[idx],0)+2.5*np.degrees(np.arccos(np.clip(c,0,1))))
        ch=self.C[idx[k]]; i=int(self.I[idx[k]])
        return ch,i,float(lat[k]),int(self.mode[ch][i]),self.T[idx[k]]

class Lines:
    """all lane and edge lines, for finding a line's neighbours across the road"""
    def __init__(self,lanes,edges):
        segs=[];kind=[];lid=[]
        for k,ls in ((0,lanes),(1,edges)):
            for j,a in enumerate(ls):
                for p,q in zip(a[:-1],a[1:]):
                    if np.hypot(*(q-p))>1e-3: segs.append(sg.LineString([p,q])); kind.append(k); lid.append((k,j))
        self.segs=segs; self.kind=np.array(kind); self.lid=lid; self.tree=shapely.STRtree(segs)
        self.dir=np.array([np.diff(np.array(s.coords),axis=0)[0]/max(s.length,1e-9) for s in segs])
    def across(self,p,t,me,reach=6.5):
        """lines crossing the normal at p (within reach), as (signed offset +left, kind); me = (kind, line id) to skip"""
        n=np.array([-t[1],t[0]]); cut=sg.LineString([p-n*reach,p+n*reach]); out=[]
        for k in self.tree.query(cut,predicate='intersects'):
            if self.lid[k]==me or abs(self.dir[k]@t)<0.85: continue
            x=cut.intersection(self.segs[k])
            if x.is_empty or x.geom_type!='Point': continue
            o=float((np.array(x.coords[0])-p)@n)
            if abs(o)>0.3: out.append((o,int(self.kind[k])))
        return out

def roles(lanes,edges,traffic,arrows,step=4.0):
    """per line, per vertex: lanes -> 'dash' or 'dyel'; edges -> 'white' or 'yellow' (lines densified to step)"""
    LN=Lines(lanes,edges); At=cKDTree(arrows[:,:2]) if len(arrows) else None
    def side_arrows(p,t,n):
        if At is None: return None,None
        res={1:[],-1:[]}
        for k in At.query_ball_point(p,50.0):
            a=arrows[k]; dv=a[:2]-p; along=dv@t; lat=dv@n; dot=a[2:4]@t
            if abs(dot)<0.85 or abs(along)>50 or not 0.6<abs(lat)<4.5: continue
            res[int(np.sign(lat))].append(np.sign(dot))
        f=lambda v: (1 if np.mean(v)>0 else -1) if v else None
        return f(res[1]),f(res[-1])
    out_l=[];out_e=[]
    for kind,ls,out in ((0,lanes,out_l),(1,edges,out_e)):
        for j,a in enumerate(ls):
            t,nv=_frame(a); lab=[]
            for i,p in enumerate(a):
                ow=traffic.owner(p,t[i])
                if kind==0:
                    al,ar=side_arrows(p,t[i],nv[i])
                    if al is not None and ar is not None: lab.append('dyel' if al!=ar else 'dash'); continue
                    if ow is None: lab.append('dash'); continue
                    ch,si,lat,md,tc=ow
                    lab.append('dyel' if md==0 and abs(lat)<1.2 else 'dash'); continue
                # an edge line: which side are its lanes on, and which way do they run?
                nb=LN.across(p,t[i],(1,j))
                near={s:[(abs(o),k) for o,k in nb if np.sign(o)==s and LANE[0]<=abs(o)<=LANE[1]] for s in (1,-1)}
                lane_side=None
                for s in (1,-1):
                    if any(k==0 for _,k in near[s]) and not any(k==0 for _,k in near[-s]): lane_side=s
                if lane_side is None and ow is not None:
                    ch,si,lat,md,tc=ow
                    if near[1] and not near[-1]: lane_side=1
                    elif near[-1] and not near[1]: lane_side=-1
                    elif abs(lat)>1.0: lane_side=-int(np.sign(lat))*int(np.sign(nv[i]@np.array([-tc[1],tc[0]])) or 1)   # towards the alignment
                if lane_side is None or ow is None: lab.append('white'); continue
                ch,si,lat,md,tc=ow
                # travel direction of those lanes, in this line's frame
                if md!=0: d=md*np.sign(tc@t[i])
                else:
                    ql=lat+lane_side*1.8*np.sign(nv[i]@np.array([-tc[1],tc[0]]))   # lane centre's offset from the alignment
                    d=(1 if ql<0 else -1)*np.sign(tc@t[i])
                # lanes to the right of the line (looking the way they run) make it their left edge
                lab.append('yellow' if lane_side*d<0 else 'white')
            out.append(np.array(lab))
    return out_l,out_e

def centre_lines(cls,chains,hw,traffic,lines,pave,step=5.0):
    """A two-way road the plans draw with no line down its middle (the alignment is drawn there instead): a double yellow
    centre line along the alignment where its pavement runs with edge lines on both sides and nothing else within
    2 m. Plan alignments count when their arrows (or lack of any) make them two-way; Overture street centrelines only
    where the plan arrows show traffic both ways, since a divided highway's carriageways are drawn one per line."""
    out=[]; shapely.prepare(pave)
    for c,v in cls.items():
        if c.startswith(('e_','y')) or hw(c,chains)>9.0 or chains[c].get('cls') in ('motorway','trunk'): continue
        P=v['P']; t,nv=_frame(P); md=traffic.mode[c]; ok=np.zeros(len(P),bool)
        isx=c.startswith('x')
        for i in range(len(P)):
            if md[i]!=0 or (isx and not traffic.both[c][i]): continue
            if not pave.contains(sg.Point(*P[i])): continue
            nb=lines.across(P[i],t[i],None,reach=hw(c,chains)+3.0)
            if any(abs(o)<2.0 for o,k in nb): continue
            L=[o for o,k in nb if k==1 and o>=2.7]; R=[o for o,k in nb if k==1 and o<=-2.7]
            ok[i]=bool(L) and bool(R)
        idx=np.where(ok)[0]
        for r in (np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(idx) else []):
            if len(r)*step>=20: out.append(P[r])
    return out

def smooth(lab,k=2):
    """majority over 2k+1 vertices, so a role does not flicker along a line"""
    lab=np.asarray(lab); n=len(lab)
    if n<3: return lab
    u=np.unique(lab); out=lab.copy()
    for i in range(n):
        w=lab[max(0,i-k):i+k+1]; c=[(w==x).sum() for x in u]; out[i]=u[int(np.argmax(c))]
    return out

def chain_lines(lines,tol=0.4,ang=0.9):
    """join line pieces that meet end to end (the plans break lines at sheet edges and curve points, in either
    direction), so dashes run evenly along them"""
    L=[np.asarray(a,float) for a in lines if len(a)>=2 and np.hypot(*np.diff(a,axis=0).T).sum()>0.5]
    if not L: return []
    E=np.array([a[0] for a in L]+[a[-1] for a in L]); tree=cKDTree(E); n=len(L); used=[False]*n; out=[]
    def head(a): v=a[-1]-a[-2]; return v/(np.hypot(*v)+1e-9)
    def grow(cur):
        while True:
            hit=None
            for e in tree.query_ball_point(cur[-1],tol):
                j=e%n
                if used[j]: continue
                b=L[j] if e<n else L[j][::-1]
                v=b[1]-b[0]; v=v/(np.hypot(*v)+1e-9)
                if head(cur)@v>ang: hit=(j,b); break
            if hit is None: return cur
            used[hit[0]]=True; cur=np.vstack([cur,hit[1][1:]])
    for i in range(n):
        if used[i]: continue
        used[i]=True; cur=grow(L[i]); cur=grow(cur[::-1])[::-1]; out.append(cur)
    return out
