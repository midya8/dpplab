import json,random,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dpplab import audit
r=Path(__file__).resolve().parents[1]
def spec(a,inc=[],exc=[],k=None):return {'L':[[str(v) for v in row] for row in a],'include':inc,'exclude':exc,'k':k}
cases=[('Identity is independent, not fixed size',spec([[1,0],[0,1]])),('Duplicate items cannot co-occur',spec([[4,4,0],[4,4,0],[0,0,1]])),('Include one duplicate',spec([[4,4,0],[4,4,0],[0,0,1]],[0])),('Impossible duplicate pair',spec([[4,4,0],[4,4,0],[0,0,1]],[0,1])),('Negative similarity still repels',spec([[2,-1],[-1,2]])),('Require exactly two items',spec([[2,1,0],[1,2,0],[0,0,1]],k=2)),('Rank-one fixed-size event impossible',spec([[1,1],[1,1]],k=2)),('Zero kernel: only empty',spec([[0,0],[0,0]])),('Exclude then include',spec([[3,1,1],[1,2,0],[1,0,2]],[1],[2])),('All items excluded',spec([[2,1],[1,2]],[],[0,1])),('All items included',spec([[2,1],[1,2]],[0,1])),('Zero item inclusion impossible',spec([[0,0],[0,2]],[0])),('Fractional kernel',spec([['1/2','1/3'],['1/3','2/3']])),('Empty-size conditioning',spec([[2,1],[1,2]],k=0)),('Fixed size after inclusion',spec([[3,1,1],[1,2,0],[1,0,2]],[0],[],2))]
rng=random.Random(20261109)
for i in range(9):
 n=3+i%3;b=[[rng.randint(-2,2) for _ in range(3)] for _ in range(n)];a=[[sum(x*y for x,y in zip(b[i],b[j])) for j in range(n)] for i in range(n)];cases.append(('Seeded Gram '+str(i),spec(a,[0] if i%3==1 else [],[n-1] if i%3==2 else [],i%3 if i%2 else None)))
out={'cases':[{'name':name,'report':audit(s)} for name,s in cases]};(r/'docs/audits.json').write_text(json.dumps(out,separators=(',',':'))+'\n');(r/'examples/duplicates.json').write_text(json.dumps(cases[1][1],indent=2)+'\n');print('Retained',len(cases),'exact reports')
