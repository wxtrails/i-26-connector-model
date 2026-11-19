"""Hand-placed fixes to pavement levels (like RIDE and BRIDGE_REACH in profiles2.py), used by meshz.py.

meshz.py gives each ground pavement vertex the level of the road whose carriageway it lies in (else the nearest edge).
Where two roads at different levels meet at the edge of a bridge, that rule can hand pavement to the wrong one, and the
triangles reaching from it to the rest of its road come out too steep: they are dropped, leaving a hole down to the
ground, or drawn as a wall slanting into the lanes. Each site gives:
- regions: (polygon, chain) pairs; ground pavement inside takes that chain's level;
- deck edges: within 1.5 m of one, ground pavement is under the bridge only if it is more than 0.3 m inside the deck
  itself, not in the metre round it that meshz.py otherwise counts as under.
"""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg

def _deck_at(plan,p):
    pt=sg.Point(*p)
    for d in plan['decks']:
        if d['poly'].contains(pt): return max(getattr(d['poly'],'geoms',[d['poly']]),key=lambda g:g.area if g.contains(pt) else 0)
    raise ValueError('no deck at %s'%(p,))

def north_end(plan,cls):
    """North end of the river bridges, where the I-240 East ramp (b2) runs down the west side of the I-26 south-bound
    bridge and dives under it at a skew.
    - The ramp's pavement reaches east to a gap that parts it from the south-bound lanes (5 m higher) and, past the gap's
      tip, to the bridge's west edge. I-26 is 44 m wide, so its carriageway reaches over that shoulder and took the
      pavement along the gap and the bridge edge to its level, leaving a wedge of dropped triangles (and a slanting wall)
      across the ramp's lanes. East of the ramp's centre line, south of a line through the gap's tip and the bridge's
      north-west corner, the pavement is the ramp's.
    - East of the bridge the north-bound lanes (whose own bridge ends 110 m further south) are on the ground at the
      deck's level, and their pavement stops at the deck's east edge. Counted as under the bridge, that edge took the
      ramp's level 7 m below, leaving a strip of holes along the bridge: the deck's east edge, and its north edge where
      the north-bound lanes go onto it, are strict deck edges."""
    D=_deck_at(plan,(285700,211370)).simplify(0.25)
    tip=np.array([285665.8,211416.2]); nw=np.array([285671.8,211414.5])        # gap's tip, bridge's north-west corner
    v=(nw-tip)/np.hypot(*(nw-tip)); left=np.array([-v[1],v[0]])
    south=sg.Polygon([tip-v*40,nw+v*40,nw+v*40-left*90,tip-v*40-left*90])
    east=sg.LineString(cls['b2']['P']).buffer(-16.0,single_sided=True)        # b2 runs north: its right is east
    side=east.intersection(south).intersection(sg.box(285600,211330,285760,211440))
    ring=sg.LineString(D.exterior.coords); a,b=ring.project(sg.Point(285688.4,211422.7)),ring.project(sg.Point(285768.0,211316.0))
    import shapely.ops as so
    edge=so.substring(ring,min(a,b),max(a,b))
    if edge.length>2*abs(a-b)+1: edge=so.substring(ring,max(a,b),ring.length)       # (the short way round)
    return [(side,'b2')],[edge]

SITES=[north_end]
def build(plan,cls):
    """(regions [(polygon, chain)], strict deck edges [LineString]) from all sites"""
    R=[];E=[]
    for f in SITES:
        r,e=f(plan,cls); R+=r; E+=e
    return R,E
