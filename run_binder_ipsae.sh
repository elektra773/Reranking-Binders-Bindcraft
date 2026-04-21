#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_binder_ipsae.sh path/to/run.conf

The config file is a bash file with KEY=VALUE assignments.
See example_ec1_run.conf for a template.
EOF
}

if [[ $# -ne 1 ]]; then
  usage
  exit 1
fi

CONFIG_PATH="$1"
if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"

maybe_init_modules() {
  if command -v module >/dev/null 2>&1; then
    return 0
  fi
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    # shellcheck disable=SC1091
    source /usr/share/Modules/init/bash
  fi
}

# shellcheck disable=SC1090
source "$CONFIG_PATH"

MODE="${MODE:-predict-and-score}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-$PROJECT_ROOT/binder_ipsae_pipeline.py}"
BOLTZ_MODULE_NAME="${BOLTZ_MODULE_NAME:-}"
TARGET_PDB="${TARGET_PDB:-}"
TARGET_PDB_CHAIN="${TARGET_PDB_CHAIN:-A}"
TARGET_CHAIN_ID="${TARGET_CHAIN_ID:-A}"
TARGET_NAME="${TARGET_NAME:-target}"
TARGET_CALCIUM_IONS="${TARGET_CALCIUM_IONS:-0}"
TARGET_MSA="${TARGET_MSA:-}"
BINDER_CSV="${BINDER_CSV:-}"
BINDER_FASTA="${BINDER_FASTA:-}"
BINDER_SEQUENCE="${BINDER_SEQUENCE:-}"
BINDER_NAME="${BINDER_NAME:-binder_1}"
BINDER_NAME_COLUMN="${BINDER_NAME_COLUMN:-Design}"
BINDER_SEQUENCE_COLUMN="${BINDER_SEQUENCE_COLUMN:-Sequence}"
BINDER_RANK_COLUMN="${BINDER_RANK_COLUMN:-Rank}"
BINDER_CHAIN_ID="${BINDER_CHAIN_ID:-B}"
BINDER_MSA="${BINDER_MSA:-}"
WORK_DIR="${WORK_DIR:-}"
SUMMARY_CSV="${SUMMARY_CSV:-}"
USE_MSA_SERVER="${USE_MSA_SERVER:-1}"
USE_POTENTIALS="${USE_POTENTIALS:-0}"
FORCE_EMPTY_MSA="${FORCE_EMPTY_MSA:-0}"
OVERRIDE="${OVERRIDE:-0}"
DEVICES="${DEVICES:-1}"
RECYCLING_STEPS="${RECYCLING_STEPS:-}"
DIFFUSION_SAMPLES="${DIFFUSION_SAMPLES:-}"
PAE_CUTOFF="${PAE_CUTOFF:-15}"
DIST_CUTOFF="${DIST_CUTOFF:-15}"
LOG_FILE="${LOG_FILE:-}"

if [[ -n "$BOLTZ_MODULE_NAME" ]]; then
  maybe_init_modules
  if ! command -v module >/dev/null 2>&1; then
    echo "Could not initialize environment modules, but BOLTZ_MODULE_NAME is set." >&2
    exit 1
  fi
  module load "$BOLTZ_MODULE_NAME"
fi

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
  echo "Pipeline script not found: $PIPELINE_SCRIPT" >&2
  exit 1
fi

if [[ "$MODE" != "score-batch" && -z "$TARGET_PDB" ]]; then
  echo "TARGET_PDB is required for MODE=$MODE" >&2
  exit 1
fi

if [[ "$MODE" != "score-batch" && -z "$WORK_DIR" ]]; then
  echo "WORK_DIR is required for MODE=$MODE" >&2
  exit 1
fi

if [[ "$MODE" != "score-batch" && -z "$BINDER_CSV" && -z "$BINDER_FASTA" && -z "$BINDER_SEQUENCE" ]]; then
  echo "Provide one of BINDER_CSV, BINDER_FASTA, or BINDER_SEQUENCE in the config." >&2
  exit 1
fi

if [[ "$MODE" != "score-batch" ]]; then
  mkdir -p "$WORK_DIR"
fi

if [[ -z "$SUMMARY_CSV" ]]; then
  if [[ "$MODE" == "score-batch" ]]; then
    echo "SUMMARY_CSV is required for MODE=score-batch" >&2
    exit 1
  fi
  SUMMARY_CSV="$WORK_DIR/ipsae_summary.csv"
fi

if [[ -z "$LOG_FILE" ]]; then
  timestamp="$(date +%Y%m%d_%H%M%S)"
  if [[ "$MODE" == "score-batch" ]]; then
    LOG_FILE="$(dirname "$SUMMARY_CSV")/run_${timestamp}.log"
  else
    LOG_FILE="$WORK_DIR/run_${timestamp}.log"
  fi
fi

mkdir -p "$(dirname "$SUMMARY_CSV")"
mkdir -p "$(dirname "$LOG_FILE")"

args=("$PIPELINE_SCRIPT" "$MODE")

if [[ "$MODE" != "score-batch" ]]; then
  args+=(
    --target-pdb "$TARGET_PDB"
    --target-pdb-chain "$TARGET_PDB_CHAIN"
    --target-chain-id "$TARGET_CHAIN_ID"
    --target-name "$TARGET_NAME"
    --target-calcium-ions "$TARGET_CALCIUM_IONS"
    --binder-chain-id "$BINDER_CHAIN_ID"
  )

  if [[ -n "$TARGET_MSA" ]]; then
    args+=(--target-msa "$TARGET_MSA")
  fi
  if [[ -n "$BINDER_MSA" ]]; then
    args+=(--binder-msa "$BINDER_MSA")
  fi
  if [[ -n "$BINDER_CSV" ]]; then
    args+=(
      --binder-csv "$BINDER_CSV"
      --binder-name-column "$BINDER_NAME_COLUMN"
      --binder-sequence-column "$BINDER_SEQUENCE_COLUMN"
      --binder-rank-column "$BINDER_RANK_COLUMN"
    )
  fi
  if [[ -n "$BINDER_FASTA" ]]; then
    args+=(--binder-fasta "$BINDER_FASTA")
  fi
  if [[ -n "$BINDER_SEQUENCE" ]]; then
    args+=(--binder-sequence "$BINDER_SEQUENCE" --binder-name "$BINDER_NAME")
  fi
  args+=(--work-dir "$WORK_DIR")
fi

if [[ "$MODE" == "predict-and-score" ]]; then
  args+=(--summary-csv "$SUMMARY_CSV" --devices "$DEVICES" --pae-cutoff "$PAE_CUTOFF" --dist-cutoff "$DIST_CUTOFF")
  if [[ "$USE_MSA_SERVER" == "1" ]]; then
    args+=(--use-msa-server)
  fi
  if [[ "$USE_POTENTIALS" == "1" ]]; then
    args+=(--use-potentials)
  fi
  if [[ "$FORCE_EMPTY_MSA" == "1" ]]; then
    args+=(--force-empty-msa)
  fi
  if [[ "$OVERRIDE" == "1" ]]; then
    args+=(--override)
  fi
  if [[ -n "$RECYCLING_STEPS" ]]; then
    args+=(--recycling-steps "$RECYCLING_STEPS")
  fi
  if [[ -n "$DIFFUSION_SAMPLES" ]]; then
    args+=(--diffusion-samples "$DIFFUSION_SAMPLES")
  fi
elif [[ "$MODE" == "prepare" ]]; then
  if [[ "$FORCE_EMPTY_MSA" == "1" ]]; then
    args+=(--force-empty-msa)
  fi
elif [[ "$MODE" == "score-batch" ]]; then
  PREDICTIONS_ROOT="${PREDICTIONS_ROOT:-}"
  if [[ -z "$PREDICTIONS_ROOT" ]]; then
    echo "PREDICTIONS_ROOT is required for MODE=score-batch" >&2
    exit 1
  fi
  args+=(--predictions-root "$PREDICTIONS_ROOT" --summary-csv "$SUMMARY_CSV" --pae-cutoff "$PAE_CUTOFF" --dist-cutoff "$DIST_CUTOFF")
fi

echo "Logging to $LOG_FILE"
echo "Running command:"
printf '  %q' "$PYTHON_BIN" "${args[@]}"
printf '\n'

"$PYTHON_BIN" "${args[@]}" 2>&1 | tee "$LOG_FILE"
