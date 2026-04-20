# Generalized Binder iPSAE Pipeline

This folder recreates the useful part of Adaptyv Bio's Nipah iPSAE notebook, but makes the target configurable so you can score your own binder designs against a different protein.

## What is included

- `binder_ipsae_pipeline.py`
  - prepares Boltz YAML inputs for one target and many binders
  - optionally runs Boltz predictions
  - scores each prediction with `ipsae.py`
  - writes a flat `ipsae_summary.csv`
- `ipsae.py`
  - vendored from the Adaptyv/Dunbrack workflow and kept as a separate script for reproducibility
- `requirements.txt`
  - minimal local dependency needed for `ipsae.py`

## Install

`ipsae.py` needs `numpy`.

```bash
python3 -m pip install -r requirements.txt
```

If you want this script to also run predictions, install Boltz separately. The current Boltz docs describe:

```bash
python3 -m pip install boltz
```

## Typical workflow

### 1. Prepare inputs from a target PDB and a FASTA of binders

This is useful when you already know the target structure or at least have a target chain sequence in a PDB.

```bash
python3 binder_ipsae_pipeline.py prepare \
  --target-pdb EC1-2/EC1.pdb \
  --target-pdb-chain A \
  --target-name CDH23_EC1 \
  --binder-fasta binders.fasta \
  --work-dir binder_ipsae_runs
```

This creates:

- `binder_ipsae_runs/inputs/*.yaml`
- `binder_ipsae_runs/binder_manifest.csv`

### 1b. Prepare inputs directly from BindCraft `final_design_stats.csv`

If your sequences already live in a BindCraft-style CSV with `Rank`, `Design`, and `Sequence` columns:

```bash
python3 binder_ipsae_pipeline.py prepare \
  --target-pdb EC1-2/EC1.pdb \
  --target-pdb-chain A \
  --target-name CDH23_EC1 \
  --binder-csv EC1-2/final_design_stats.csv \
  --work-dir binder_ipsae_runs
```

This keeps the ranking order from the CSV and names jobs like `001_1_CDH23EC1_...`.

### 2. Run Boltz and score in one step

If you want Boltz to build the complexes and then score them immediately:

```bash
python3 binder_ipsae_pipeline.py predict-and-score \
  --target-pdb EC1-2/EC1.pdb \
  --target-pdb-chain A \
  --target-name CDH23_EC1 \
  --binder-fasta binders.fasta \
  --work-dir binder_ipsae_runs \
  --use-msa-server \
  --summary-csv binder_ipsae_runs/ipsae_summary.csv
```

Optional flags you may also want:

- `--target-msa path/to/target.a3m`
- `--use-potentials`
- `--override`
- `--recycling-steps 10`
- `--diffusion-samples 5`

### 3. Score predictions you already have

For a single prediction:

```bash
python3 binder_ipsae_pipeline.py score \
  --pae path/to/pae_example_model_0.npz \
  --structure path/to/example_model_0.cif \
  --output-json score.json \
  --output-csv score.csv
```

For an entire Boltz output folder:

```bash
python3 binder_ipsae_pipeline.py score-batch \
  --predictions-root binder_ipsae_runs/boltz_output \
  --summary-csv binder_ipsae_runs/ipsae_summary.csv
```

## Input assumptions

- This wrapper is built for the common case of one target protein chain and one binder chain.
- The target is placed on chain `A` and the binder on chain `B` by default.
- If you do not provide an MSA path, the YAML uses `msa: empty`. In that case you should usually run Boltz with `--use-msa-server`.

## Output summary

The flat CSV contains one `max` row and one `min` row per chain pair.

- `max` matches the competition-style maximum inter-chain iPSAE summary.
- `min` is the minimum of the two asymmetric directions, which is useful as a stricter filter.

## Provenance

This setup is adapted from:

- Adaptyv Bio's Nipah workflow: [adaptyvbio/nipah_ipsae_pipeline](https://github.com/adaptyvbio/nipah_ipsae_pipeline)
- The bundled scorer script: [DunbrackLab/IPSAE](https://github.com/DunbrackLab/IPSAE)
- Boltz prediction docs: [jwohlwend/boltz prediction docs](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)
