import json,sys,numpy as np,collections
f=sys.argv[1]; R=[json.loads(l) for l in open(f)]
by=collections.defaultdict(list)
for r in R: by[r['sample']].append(r)
out={}
for s,L in by.items():
    rv=[r['reverse'] for r in L if 'reverse' in r]
    a=lambda k:np.array([x[k] for x in rv],float)
    out[s]=dict(n=len(L),n_err=len(L)-len(rv),n_pass=int(sum(x['pair_pass'] for x in rv)),
      n_in_files=int(sum(x['in_files'] for x in rv)),
      frac_on_median=float(np.median(a('frac_on'))),frac_on_min=float(a('frac_on').min()),n_frac_ok=int((a('frac_on')>=.95).sum()),
      median_of_median_normal=float(np.median(a('median_normal_coarse_vx'))),
      p95_normal_median=float(np.median(a('p95_normal_coarse_vx'))),
      inplane_p95_median=float(np.median(a('inplane_p95_px'))),n_inplane_ok=int((a('inplane_p95_px')<=.003).sum()))
print(json.dumps(out,indent=1))
