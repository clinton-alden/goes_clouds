#!/usr/bin/env bash
set -euo pipefail

NB=${1:-analysis/11_compare_sw_tsi.ipynb}
OUTDIR=${2:-}

python - <<'PY'
import importlib.util
import subprocess
import sys

if importlib.util.find_spec('pyproj') is None:
    print('pyproj not found; installing with pip for this environment...')
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'pyproj'])
PY

export MPLCONFIGDIR=${MPLCONFIGDIR:-/glade/derecho/scratch/cdalden/tmp/mplconfig}
mkdir -p "$MPLCONFIGDIR"

if command -v jupyter-nbconvert >/dev/null 2>&1; then
    if jupyter-nbconvert --to notebook --execute "$NB" --inplace ${OUTDIR:+--output-dir "$OUTDIR"}; then
        echo "nbconvert execution complete: $NB"
        exit 0
    fi
fi

python - <<PY
import importlib.util, json, traceback
from pathlib import Path

if importlib.util.find_spec('nbconvert'):
    import subprocess, shlex, sys
    nb_path = Path('$NB')
    cmd = ['python', '-m', 'nbconvert', '--to', 'notebook', '--execute', str(nb_path), '--inplace']
    if '$OUTDIR':
        cmd.extend(['--output-dir', '$OUTDIR'])
    rc = subprocess.call(cmd)
    if rc == 0:
        print(f'nbconvert execution complete: {nb_path}')
        raise SystemExit(0)

# Fallback: execute notebook cells as plain script
nb_path = Path('$NB')
nb = json.loads(nb_path.read_text())
code = '\n\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code')

ns = {'__name__': '__main__'}
try:
    exec(compile(code, str(nb_path), 'exec'), ns)
    print(f'Execution completed via fallback script: {nb_path}')
except Exception:
    traceback.print_exc()
    raise SystemExit(1)
PY
