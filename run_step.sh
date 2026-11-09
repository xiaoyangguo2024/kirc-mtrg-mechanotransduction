#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/Rscript" ]; then
  RSCRIPT="${CONDA_PREFIX}/bin/Rscript"
elif [ -x "../BLCA_Rloop/conda_env/bin/Rscript" ]; then
  export PATH="$(cd ../BLCA_Rloop/conda_env/bin && pwd):$PATH"
  RSCRIPT="Rscript"
else
  RSCRIPT="Rscript"
fi
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
mkdir -p logs
STEP="$1"
echo "[$(date '+%F %T')] Starting ${STEP}"
"${RSCRIPT}" "R/${STEP}.R" > "logs/${STEP}.log" 2>&1
echo "[$(date '+%F %T')] Finished ${STEP} (exit $?)"
