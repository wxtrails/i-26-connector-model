import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from mosaic import *
import json
SX0,SY0,SX1,SY1=279480.0,204480.0,288600.0,214080.0
def adj(a):
    f=np.moveaxis(a,0,-1).astype(np.float32)/255; nod=(a.sum(0)==0)
    g=np.clip((f-0.02)/0.72,0,1)**0.78
    g[nod]=(0.36,0.42,0.30)
    return (g*255).astype(np.uint8)
def base():
    a=mosaic(SX0,SY0,SX1,SY1,2.4); Image.fromarray(adj(a)).save(W+'/scene/base.jpg',quality=80,optimize=True); print('base',a.shape)
def atlas():
    T=np.load(W+'/terrain_tiles.npz'); tiles=T['tiles']; TILE=240.0; RES=0.6; N=400; GUT=4; CELL=N+2*GUT; PER=10
    blocks={}
    for k,(i,j) in enumerate(tiles): blocks.setdefault((i//4,j//4),[]).append((k,i,j))
    nat=(len(tiles)+PER*PER-1)//(PER*PER)
    at=[np.zeros((CELL*PER,CELL*PER,3),np.uint8) for _ in range(nat)]
    for (bi,bj),lst in blocks.items():
        x0=SX0+bi*4*TILE-GUT*RES; y0=SY0+bj*4*TILE-GUT*RES
        a=adj(mosaic(x0,y0,x0+4*TILE+2*GUT*RES,y0+4*TILE+2*GUT*RES,RES))   # rows from top (y max)
        H=a.shape[0]
        for k,i,j in lst:
            ci=(i-bi*4)*N; cj=(j-bj*4)*N     # px offset from left / from bottom
            cell=a[H-(cj+CELL):H-cj, ci:ci+CELL]
            A=at[k//(PER*PER)]; q=k%(PER*PER); r,c=q//PER,q%PER
            A[r*CELL:(r+1)*CELL,c*CELL:(c+1)*CELL]=cell
        print('block',bi,bj,len(lst),flush=True)
    for n,A in enumerate(at): Image.fromarray(A).save(W+f'/scene/atlas{n}.jpg',quality=76,optimize=True)
    # record which terrain tile each cell belongs to, so export.py can tell when the atlas no longer matches terrain_tiles.npz
    json.dump({'n':nat,'per':PER,'cell':CELL,'gut':GUT,'tile':N,'tiles':tiles.tolist()},open(W+'/scene/atlas.json','w'))
if __name__=='__main__':
    what=sys.argv[1:] or ['base','atlas']      # e.g. `python imagery.py atlas` re-cuts only the atlases
    if 'base' in what: base()
    if 'atlas' in what: atlas()
    print('IMAGERY DONE')
