import argparse,json,sys
from pathlib import Path
from . import audit

def strict(pairs):
    out={}
    for k,v in pairs:
        if k in out:raise ValueError('Duplicate JSON key')
        out[k]=v
    return out

def main():
    p=argparse.ArgumentParser(description='Exact finite L-ensemble audit');p.add_argument('--spec',required=True);p.add_argument('--json',required=True);args=p.parse_args()
    try:
        data=Path(args.spec).read_bytes()
        if len(data)>32768:raise ValueError('Input limit 32KiB')
        spec=json.loads(data,object_pairs_hook=strict,parse_constant=lambda _:(_ for _ in ()).throw(ValueError('Nonfinite JSON token')))
        result=audit(spec);Path(args.json).write_text(json.dumps(result,indent=2)+'\n');print('possible' if result['selected_mean'] is not None else 'impossible conditioning event');return 0
    except (ValueError,OSError) as e:print('ERROR '+str(e),file=sys.stderr);return 2
