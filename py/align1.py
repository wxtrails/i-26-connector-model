import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely.geometry as sg, shapely.ops as so
allL=[]
for n in (1,2,3):
    L=layer_lines(n,'Prop Horizontal Alignment',('s',))
    L=[l for l in L if len(l)>=2]
    lens=[np.hypot(*np.diff(l,axis=0).T).sum() for l in L]
    print(n,len(L),'total len',sum(lens),'long',sorted([round(x) for x in lens if x>60],reverse=True)[:40])
    for l,ln in zip(L,lens): allL.append((n,l,ln))
pickle.dump(allL,open(W+'/align_raw.pkl','wb'))
