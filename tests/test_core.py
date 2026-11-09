import copy,itertools,json,random,subprocess,sys,tempfile,unittest
from fractions import Fraction as F
from pathlib import Path
from dpplab import audit
from dpplab.core import determinant,inverse

def refdet(a):
    total=F(0)
    for p in itertools.permutations(range(len(a))):
        x=F((-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p))))
        for i,j in enumerate(p):x*=F(a[i][j])
        total+=x
    return total

def principal(a,s):return [[a[i][j] for j in s] for i in s]
def subsets(n):return [list(s) for k in range(n+1) for s in itertools.combinations(range(n),k)]
def gram(rng,n):
    b=[[rng.randint(-3,3) for _ in range(rng.randint(1,1)+2)] for _ in range(n)]
    return [[sum(x*y for x,y in zip(b[i],b[j])) for j in range(n)] for i in range(n)]
def spec(a,inc=[],exc=[],k=None):return {'L':[[str(v) for v in r] for r in a],'include':inc,'exclude':exc,'k':k}

class ExactTests(unittest.TestCase):
    def test_independent_determinants(self):
        rng=random.Random(71)
        for n in range(6):
            for _ in range(20):
                a=[[F(rng.randint(-4,4),rng.randint(1,4)) for _ in range(n)] for _ in range(n)]
                self.assertEqual(determinant(a),refdet(a))
    def test_integral_json_numbers_and_fixed_k_kernel_distinction(self):
        a=spec([[1,0],[0,1]],inc=[0.0],k=1.0);r=audit(a)
        self.assertEqual(r['selected_inclusions'],['1','0'])
        b=audit(spec([[1,0],[0,1]],k=1))
        self.assertEqual(b['selected_inclusions'],['1/2','1/2'])
        self.assertEqual(b['rows'][-1]['conditional_probability'],'0')
        self.assertEqual(b['rows'][-1]['inclusion_determinant'],'1/4')
        self.assertEqual(b['event_conditional_K'],[['1/2','0'],['0','1/2']])
    def test_inverse_multiplication(self):
        rng=random.Random(72)
        for _ in range(50):
            n=rng.randint(1,5);a=gram(rng,n);a=[[F(a[i][j]+(i==j)) for j in range(n)] for i in range(n)];b=inverse(a)
            for i in range(n):
                for j in range(n):self.assertEqual(sum(a[i][k]*b[k][j] for k in range(n)),F(i==j))
    def test_140_gram_all_minor_inclusion_size_oracles(self):
        rng=random.Random(73)
        for _ in range(140):
            n=rng.randint(1,5);a=gram(rng,n);r=audit(spec(a));ss=subsets(n);weights=[refdet(principal(a,s)) for s in ss];z=sum(weights);K=[[F(x) for x in row] for row in r['K']]
            self.assertEqual(F(r['normalizer']),z)
            for row,s,w in zip(r['rows'],ss,weights):
                self.assertEqual(row['subset'],s);self.assertEqual(F(row['determinant']),w);self.assertEqual(F(row['probability']),w/z);p=sum(v for t,v in zip(ss,weights) if set(s)<=set(t))/z
                self.assertEqual(F(row['inclusion_sum']),p);self.assertEqual(refdet(principal(K,s)),p)
            mean=sum(F(len(s))*w/z for s,w in zip(ss,weights));var=sum(F(len(s)**2)*w/z for s,w in zip(ss,weights))-mean*mean
            self.assertEqual(F(r['baseline_mean']),mean);self.assertEqual(F(r['baseline_variance']),var)
            self.assertEqual([F(v) for v in r['selected_cardinality']],[sum((w/z for s,w in zip(ss,weights) if len(s)==k),F(0)) for k in range(n+1)])
    def test_100_conditional_oracles_and_schur(self):
        rng=random.Random(74)
        for _ in range(100):
            n=rng.randint(1,5);a=gram(rng,n);state=[rng.randrange(3) for _ in range(n)];inc=[i for i,x in enumerate(state) if x==1];exc=[i for i,x in enumerate(state) if x==2];k=rng.choice([None,*range(n+1)]);r=audit(spec(a,inc,exc,k));ss=subsets(n);w=[refdet(principal(a,s)) for s in ss];event=[set(inc)<=set(s) and not set(exc)&set(s) for s in ss];match=[e and (k is None or len(s)==k) for e,s in zip(event,ss)];ez=sum(v for v,e in zip(w,event) if e);sz=sum(v for v,e in zip(w,match) if e)
            self.assertEqual(F(r['event_weight']),ez);self.assertEqual(F(r['selected_weight']),sz)
            self.assertEqual([row['conditional_probability'] for row in r['rows']],[str(v/sz if m else F(0)) if sz else None for v,m in zip(w,match)])
            if not ez:self.assertIsNone(r['event_conditional_L']);continue
            rem=r['remaining'];c=[[F(x) for x in row] for row in r['event_conditional_L']];den=refdet([[c[i][j]+(i==j) for j in range(len(c))] for i in range(len(c))]);ck=[[F(x) for x in row] for row in r['event_conditional_K']]
            for t in subsets(len(rem)):
                s=sorted(inc+[rem[i] for i in t]);idx=ss.index(s)
                self.assertEqual(refdet(principal(c,t))/den,w[idx]/ez)
                p=sum(v/ez for s,v,e in zip(ss,w,event) if e and all(rem[i] in s for i in t))
                self.assertEqual(refdet(principal(ck,t)),p)
    def test_scale_invariance_only_fixed_size(self):
        rng=random.Random(75)
        for _ in range(40):
            n=4;a=gram(rng,n);b=[[F(7,3)*x for x in row] for row in a];k=rng.randrange(5);r=audit(spec(a,k=k));s=audit(spec(b,k=k));self.assertEqual([x['conditional_probability'] for x in r['rows']],[x['conditional_probability'] for x in s['rows']])
        self.assertNotEqual(audit(spec([[1]]))['baseline_mean'],audit(spec([[7]]))['baseline_mean'])
    def test_sign_conjugation_invariance(self):
        rng=random.Random(76)
        for _ in range(20):
            a=gram(rng,4);sg=[rng.choice([-1,1]) for _ in a];b=[[x*sg[i]*sg[j] for j,x in enumerate(row)] for i,row in enumerate(a)];r=audit(spec(a));s=audit(spec(b));self.assertEqual(r['rows'],s['rows']);self.assertEqual(r['pairs'],s['pairs'])
    def test_diagonal_independent_bernoulli(self):
        vals=[F(1,2),F(2),F(0),F(4,3)];r=audit(spec([[v if i==j else 0 for j in range(4)] for i,v in enumerate(vals)]))
        for row in r['rows']:
            p=F(1)
            for i,v in enumerate(vals):p*=v/(1+v) if i in row['subset'] else 1/(1+v)
            self.assertEqual(F(row['probability']),p)
    def test_duplicate_rank_one_and_impossible(self):
        for k in [None,0,1,2]:
            r=audit(spec([[4,4],[4,4]],k=k));self.assertEqual(r['rows'][-1]['determinant'],'0')
            if k==2:self.assertIsNone(r['selected_mean']);self.assertEqual(r['selected_modes'],[])
        self.assertIsNone(audit(spec([[4,4],[4,4]],[0,1]))['event_conditional_L'])
    def test_empty_determinant(self):self.assertEqual(determinant([]),1);self.assertEqual(inverse([]),[])
    def test_all_excluded(self):
        r=audit(spec([[2,1],[1,2]],[],[0,1]));self.assertEqual(r['selected_mean'],'0');self.assertEqual(r['event_conditional_L'],[]);self.assertEqual(r['event_conditional_K'],[])
    def test_all_included(self):
        r=audit(spec([[2,1],[1,2]],[0,1]));self.assertEqual(r['selected_mean'],'2');self.assertEqual(r['selected_inclusions'],['1','1']);self.assertEqual(r['event_conditional_L'],[])
    def test_k_is_total_not_remaining_size(self):
        r=audit(spec([[1,0],[0,1]],[0],[],1));self.assertEqual(r['selected_inclusions'],['1','0'])
    def test_modes_all_ties(self):self.assertEqual(len(audit(spec([[1,0],[0,1]]))['modes']),4)
    def test_pairs_from_distribution(self):
        r=audit(spec([[3,-1,1],[-1,2,0],[1,0,2]]));K=r['K']
        for p in r['pairs']:
            joint=sum(F(x['probability']) for x in r['rows'] if p['i'] in x['subset'] and p['j'] in x['subset']);self.assertEqual(joint,F(p['joint']));self.assertEqual(joint-F(K[p['i']][p['i']])*F(K[p['j']][p['j']]),F(p['covariance']));self.assertLessEqual(F(p['covariance']),0)
    def test_psd_not_diagonal_or_leading_only(self):
        for a in [[[1,2],[2,1]],[[1,0,0],[0,1,2],[0,2,1]],[[0,1],[1,0]]]:
            with self.assertRaisesRegex(ValueError,'semidefinite'):audit(spec(a))
    def test_asymmetric_refused(self):
        with self.assertRaisesRegex(ValueError,'symmetric'):audit(spec([[1,1],[0,1]]))
    def test_strict_rational_inputs(self):
        for x in [1,1.0,True,None,'1.0','01','1/0','1/-2','NaN','1e2','10001','1/10001','1 '*20]:
            with self.assertRaises(ValueError):audit({'L':[[x]],'include':[],'exclude':[],'k':None})
        self.assertEqual(audit(spec([['20000/2']]))['normalizer'],'10001')
    def test_index_validation(self):
        for key in ['include','exclude']:
            for val in [[True],[1.2],[-1],[2],[0,0],None]:
                s=spec([[1,0],[0,1]]);s[key]=val
                with self.assertRaises(ValueError):audit(s)
    def test_schema_and_dimensions(self):
        s=spec([[1]])
        for x in [None,[],{},dict(s,extra=0),dict(s,L=[]),dict(s,L=[['1','0']]),dict(s,L=[['1']]*6)]:
            with self.assertRaises(ValueError):audit(x)
        for k in [True,1.2,-1,2,'1']:
            with self.assertRaises(ValueError):audit(dict(s,k=k))
        with self.assertRaises(ValueError):audit(dict(s,include=[0],exclude=[0]))
    def test_n5_boundary(self):
        r=audit(spec([[F(i==j) for j in range(5)] for i in range(5)]));self.assertEqual(len(r['rows']),32);self.assertEqual(r['normalizer'],'32')
    def test_no_mutation(self):
        s=spec([[2,1],[1,2]],[1]);before=copy.deepcopy(s);audit(s);self.assertEqual(s,before)
    def test_retained_artifacts(self):
        f=json.loads((Path(__file__).resolve().parents[1]/'docs/audits.json').read_text())
        for c in f['cases']:self.assertEqual(audit(c['report']['spec']),c['report'])
    def test_cli_strict_and_impossible_success(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'in.json';out=Path(d)/'out.json';p.write_text(json.dumps(spec([[0]],[0])))
            cmd=[sys.executable,'-m','dpplab','--spec',str(p),'--json',str(out)];self.assertEqual(subprocess.run(cmd,capture_output=True).returncode,0);self.assertIsNone(json.loads(out.read_text())['selected_mean'])
            for text in ['{"L":[],"L":[]}', '{"L":NaN}', ' '*32769]:
                p.write_text(text);self.assertEqual(subprocess.run(cmd,capture_output=True).returncode,2)

if __name__=='__main__':unittest.main()
