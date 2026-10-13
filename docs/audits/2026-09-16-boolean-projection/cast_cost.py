"""Paired warm multiplication timing, separate from the RSS experiments.

Both operands are retained here to reduce allocation-history and host-drift
differences. This experiment therefore makes no memory-saving claim.
"""

import json
import random
import statistics
import time
from pathlib import Path
import benchmark as b

b.imports()
np=b.np
records=[]
for q in (1,100,1000):
    n=b.BUDGET//(8*q)
    m=8 if q==1 else 128
    for kind in ('dense','csr'):
        op=b.operator(m,n,kind,.01)
        packed,_,_,_=b.payload(n,q,.1)
        boolean=np.unpackbits(packed,axis=0,count=n,bitorder='little').view(bool)
        explicit=boolean.astype(np.float64)
        # Warm kernels and allocation pools before the paired measurements.
        op@boolean;op@explicit
        for blocks in (1,4):
            samples={'explicit64':[],'native':[]}
            cpu={'explicit64':[],'native':[]}
            rng=random.Random(883)
            for _ in range(25):
                order=['explicit64','native'];rng.shuffle(order)
                for variant in order:
                    value=explicit if variant=='explicit64' else boolean
                    wall=time.perf_counter();start=time.process_time()
                    for _ in range(blocks):
                        product=op@value
                        del product
                    cpu[variant].append(time.process_time()-start)
                    samples[variant].append(time.perf_counter()-wall)
            records.append(dict(q=q,n=n,m=m,operator=kind,blocks=blocks,
                                seconds=samples,cpu_seconds=cpu,
                                median_seconds={key:statistics.median(v) for key,v in samples.items()}))
        del op,boolean,explicit,packed
root=Path(__file__).resolve().parent
(root/'cast-cost.json').write_text(json.dumps(records,indent=2))
for r in records:
    a=r['median_seconds']['explicit64'];v=r['median_seconds']['native']
    print(r['q'],r['operator'],r['blocks'],round(a*1000,3),round(v*1000,3),round(v/a,2))
