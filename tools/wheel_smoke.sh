#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)";tmp="$(mktemp -d)";trap 'rm -rf "$tmp"' EXIT
python -m pip wheel --no-deps --wheel-dir "$tmp/dist" "$root" > "$tmp/build.log"
python -m venv "$tmp/env"
"$tmp/env/bin/python" -m pip install --no-deps "$tmp"/dist/*.whl > "$tmp/install.log"
cd "$tmp"
"$tmp/env/bin/python" -I "$root/tools/verify_installed.py" "$root"
"$tmp/env/bin/dpplab" --spec "$root/examples/duplicates.json" --json cli.json
"$tmp/env/bin/python" -I - "$root" <<'PY'
import json,sys
from pathlib import Path
from dpplab import audit
r=Path(sys.argv[1]);assert json.load(open('cli.json'))==audit(json.load(open(r/'examples/duplicates.json')))
assert json.load(open('cli.json'))['normalizer']=='18'
print('Installed console CLI, all-subset probabilities, all 24 reports and four module origins passed')
PY
