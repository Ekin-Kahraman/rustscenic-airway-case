# pySCENIC reference environment

This directory pins the environment used for the AUCell head-to-head reference in
`scripts/03_headtohead_pyscenic_aucell.py`.

The working stack is intentionally old because `pyscenic==0.12.1` and
`ctxcore==0.2.0` still depend on `pkg_resources`, older `numpy`/`numba`, and a
Python version that can also install `rustscenic>=0.4.4`.

Verified locally on 2026-05-15 with Python 3.10:

```bash
python -m venv .venv-pyscenic-reference
. .venv-pyscenic-reference/bin/activate
python -m pip install --upgrade pip
pip install -r reference/requirements-pyscenic-aucell.txt
python reference/smoke_pyscenic_reference.py
```

The smoke test checks three things:

1. `pyscenic.aucell` imports under the pinned stack.
2. `rustscenic.aucell` imports in the same environment.
3. Both AUCell implementations rank a tiny synthetic basal/ciliated matrix in
   the same direction.

This does not prove that pySCENIC's GRN stage is healthy. The case study uses
RustScenic GRN adjacencies on both sides so the comparison isolates the AUCell
kernel.
