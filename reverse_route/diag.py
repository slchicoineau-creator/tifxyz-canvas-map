"""Post-hoc descriptive diagnostic (NOT a gate): what separates off-surface coarse vertices on 0139? Uses reverse.py internals."""
import json,sys,numpy as np
from pathlib import Path
from scipy.ndimage import distance_transform_edt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from meshpair.tifxyz import read_tifxyz
from meshpair.catalog import coarse_to_fine_affine
from meshpair import reverse as R, conventions as cv
P=[p for p in json.load(open('out/pairs.json'))['pairs'] if p['sample']=='PHerc0139'][:6]
res=[]
for p in P:
    f,c=p['fine'],p['coarse']
    F=read_tifxyz('cache/meshes/'+f['key']);C=read_tifxyz('cache/meshes/'+c['key'])
    gF,gC=F.geom(f['px']),C.geom(c['px']);M=coarse_to_fine_affine(p)
    v=C.valid; rr,cc=np.nonzero(v)
    pCy=gC.grid_to_canvas(rr.astype(float),0);pCx=gC.grid_to_canvas(cc.astype(float),1)
    X=C.points[rr,cc].astype(float);Y=X@M[:,:3].T+M[:,3]
    gy0=gF.canvas_to_grid(cv.exact_canvas_map_inv(pCy,gF,gC,0),0);gx0=gF.canvas_to_grid(cv.exact_canvas_map_inv(pCx,gF,gC,1),1)
    gy,gx,nd,tg,ok=R._project(F.points,F.valid,Y,gy0,gx0)
    nd=np.abs(nd)*f['px']/c['px'];conv=ok&(tg<1e-3)
    dist=distance_transform_edt(v)[rr,cc]
    # position of the lookup inside its fine cell and inside the source-like 4x pattern
    fy=gy-np.floor(gy);fx=gx-np.floor(gx)
    off=conv&(nd>0.05)
    res.append(dict(seg=p['segment'],n=int(len(rr)),frac_off_all=float((nd>0.05).mean()),frac_notconv=float((~conv).mean()),
      off_frac_within3cells_of_edge=float((dist[off]<=3).mean()) if off.any() else None,
      all_frac_within3cells_of_edge=float((dist<=3).mean()),
      frac_off_interior_gt3=float((nd[dist>3]>0.05).mean()),
      nd_q=[float(np.quantile(nd[conv],q)) for q in (.5,.8,.9,.95,.99)]))
    print(res[-1],flush=True)
