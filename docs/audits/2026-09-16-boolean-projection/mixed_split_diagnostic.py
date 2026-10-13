"""Separate grouping-induced rounding from annotation operand precision."""

import json
from pathlib import Path
import benchmark as b

b.imports()
np=b.np
n,q,m=20971,100,128
op=b.operator(m,n,'dense',.01)
packed,continuous,bi,ci=b.payload(n,q,.1,.1)
binary=np.unpackbits(packed,axis=0,count=n,bitorder='little').view(bool)
values=np.empty((n,q),np.float32)
values[:,bi],values[:,ci]=binary,continuous
full=op@values.astype(np.float64)
native=np.empty_like(full);double=np.empty_like(full)
native[:,bi],native[:,ci]=op@binary,op@continuous
double[:,bi],double[:,ci]=op@binary.astype(np.float64),op@continuous.astype(np.float64)
result={
    'shape':[m,n,q],'continuous_fraction':.1,
    'native_split_vs_double_split_max_abs':float(np.max(np.abs(native-double))),
    'double_split_vs_full_max_abs':float(np.max(np.abs(double-full))),
    'double_split_vs_full_tolerance_violations':int(np.count_nonzero(~np.isclose(double,full,rtol=1e-14,atol=1e-14))),
    'rtol':1e-14,'atol':1e-14,
}
(Path(__file__).resolve().parent/'mixed-split-diagnostic.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
