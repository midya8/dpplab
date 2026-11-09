import importlib,json,sys
from pathlib import Path
from dpplab import audit
root=Path(sys.argv[1]).resolve();origins={}
for name in ['dpplab','dpplab.core','dpplab.cli','dpplab.__main__']:
    if name.endswith('__main__'):
        import importlib.util
        path=Path(importlib.util.find_spec(name).origin).resolve()
    else:path=Path(importlib.import_module(name).__file__).resolve()
    assert root not in path.parents,(name,path);assert 'site-packages' in str(path),(name,path);origins[name]=str(path)
f=json.loads((root/'docs/audits.json').read_text())
for c in f['cases']:assert audit(c['report']['spec'])==c['report']
print(json.dumps({'cases':len(f['cases']),'module_origins':origins}))
