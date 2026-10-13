"""Read-only local chromosome-window sizing and dtype checks for the report."""

import json
from pathlib import Path
import numpy as np
from scipy import sparse
from ldsc._kernel.ldscore import get_block_lefts

root=Path(__file__).resolve().parent
repo=root.parents[2]
bim=repo.parent/'resources/example_1kg_phase3_PLINK/1000G.EUR.QC.22.bim'
cm=np.loadtxt(bim,usecols=2)
left=get_block_lefts(cm,1).astype(np.int64)
widths=(np.ceil((np.arange(len(cm))-left)/128)*128).astype(np.int64)
sizes={}
for q in (1,54,100,1000):
    tile=16*1024**2//(8*q)
    sizes[str(q)]={
        'tile_rows':tile,
        'full_float64_tile_bytes':tile*q*8,
        'full_bool_tile_bytes':tile*q,
        'plink_forward_operand_bytes':128*q*8,
        'reverse_tiles_per_window_quantiles':np.quantile(np.ceil(widths/float(tile)),[0,.5,.9,.95,.99,1]).tolist(),
    }
promotions=[]
for optype in (np.float32,np.float64):
    for annotype in (bool,np.float32,np.float64):
        for csr in (False,True):
            op=np.eye(3,dtype=optype)
            if csr: op=sparse.csr_matrix(op)
            out=op@np.ones((3,2),dtype=annotype)
            promotions.append({'operator':np.dtype(optype).name,'annotations':np.dtype(annotype).name,'csr':csr,'result':str(out.dtype)})
result={
    'bim_path':str(bim),'bim_rows':len(cm),'cm_non_decreasing':bool(np.all(np.diff(cm)>=0)),
    'window_cm':1,'snp_batch_size':128,
    'rounded_left_window_rows_quantiles':dict(zip(('min','median','p90','p95','p99','max'),np.quantile(widths,[0,.5,.9,.95,.99,1]).tolist())),
    'warning':'BIM-only proxy before annotation/MAF/output filtering, not a full workflow profile.',
    'budgets':sizes,'dtype_promotions':promotions,
    'chr22_accumulator_1000_queries_53_baselines_weight_bytes':17380*1054*8,
    'all_output_accumulator_1000_queries_53_baselines_weight_bytes':1134471*1054*8,
}
(root/'context.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
