# Generalized Binder iPSAE Pipeline

This folder recreates the useful part of Adaptyv Bio's Nipah iPSAE notebook, but makes the target configurable so you can score your own binder designs against a different protein.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/elektra773/Reranking-Binders-Bindcraft/blob/main/Binder_iPSAE_Colab.ipynb)

## Start here

- GitHub repo: [elektra773/Reranking-Binders-Bindcraft](https://github.com/elektra773/Reranking-Binders-Bindcraft)
- Open the notebook in Colab: [Binder_iPSAE_Colab.ipynb](https://colab.research.google.com/github/elektra773/Reranking-Binders-Bindcraft/blob/main/Binder_iPSAE_Colab.ipynb)
- Local source-of-truth folder: `/Users/elektramakris/Desktop/SOTOMAYOR/Thesis/Rescoring Binders/binder-ipsae-pipeline`

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

## EC1 quick start

For your current CDH23 EC1 binder set, the easiest local entry point is the BindCraft CSV you already have:

```bash
python3 binder_ipsae_pipeline.py prepare \
  --target-pdb ../EC1-2/EC1.pdb \
  --target-pdb-chain A \
  --target-name CDH23_EC1 \
  --target-calcium-ions 3 \
  --binder-csv ../EC1-2/final_design_stats.csv \
  --work-dir binder_ipsae_runs
```

That command builds one Boltz input YAML per ranked binder and keeps the rank and design name from `final_design_stats.csv`.

## Typical workflow

### Optional: run from a shell config file

If you prefer a reusable launcher instead of a long command, use:

- `run_binder_ipsae.sh`
- `example_ec1_run.conf`

Copy the example config, edit the paths and options, then run:

```bash
cp example_ec1_run.conf my_ec1_run.conf
./run_binder_ipsae.sh my_ec1_run.conf
```

This launcher:

- can `module load boltz` for shared Linux installs
- writes a log file with `tee`
- keeps all of your main settings in one editable config file
- supports `prepare`, `predict-and-score`, and `score-batch`

### 1. Prepare inputs from a target PDB and a FASTA of binders

This is useful when you already know the target structure or at least have a target chain sequence in a PDB.

```bash
python3 binder_ipsae_pipeline.py prepare \
  --target-pdb EC1-2/EC1.pdb \
  --target-pdb-chain A \
  --target-name CDH23_EC1 \
  --target-calcium-ions 3 \
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
  --target-calcium-ions 3 \
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
  --target-calcium-ions 3 \
  --binder-csv EC1-2/final_design_stats.csv \
  --work-dir binder_ipsae_runs \
  --use-msa-server \
  --summary-csv binder_ipsae_runs/ipsae_summary.csv
```

Optional flags you may also want:

- `--target-msa path/to/target.a3m`
- `--use-potentials`
- `--force-empty-msa`
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

## Modes

- `prepare`
  - creates YAML inputs and a manifest only
  - use this when you want to inspect jobs before spending GPU time
- `predict-and-score`
  - creates inputs, runs Boltz, then runs `ipsae.py`
  - this is the full end-to-end mode
- `score`
  - scores one existing prediction from a structure file plus a matching PAE/confidence file
- `score-batch`
  - scans an existing Boltz output tree and scores every prediction it finds

In the Colab notebook the matching mode names are:

- `prepare`
- `predict_and_score`
- `score_existing_predictions`

## Input assumptions

- This wrapper is built for the common case of one target protein chain and one binder chain.
- The target is placed on chain `A` and the binder on chain `B` by default.
- If you do not provide an MSA path, the YAML now omits the `msa` field so Boltz can use `--use-msa-server` normally.
- Use `--force-empty-msa` only when you intentionally want single-sequence mode.
- If your target needs calcium ions, add `--target-calcium-ions N` and the YAML will append `N` ligand entries with CCD code `CA`.

## MSA and potentials

- `--target-msa`
  - path to a precomputed target-chain alignment file, usually `.a3m`
- `--binder-msa`
  - path to a precomputed binder-chain alignment file
- `--use-msa-server`
  - tells Boltz to build MSAs automatically when you do not already have them
  - this is the simplest option for most first runs
- `--force-empty-msa`
  - explicitly writes `msa: empty` into the YAML
  - use this only when you intentionally want single-sequence mode
- `--use-potentials`
  - passes Boltz's inference-time potentials flag
  - this can improve physical plausibility, but is usually slower than a plain run
- `--target-calcium-ions`
  - appends calcium ions as Boltz ligands with CCD code `CA`
  - useful when your target structure depends on bound calcium

For your first EC1 run, a good default is:

- leave both MSA paths blank
- turn on `--use-msa-server`
- leave `--use-potentials` off for the first pass
- optionally rerun only the best subset with `--use-potentials`

## Output summary

The flat CSV contains one `max` row and one `min` row per chain pair.

- `max` matches the competition-style maximum inter-chain iPSAE summary.
- `min` is the minimum of the two asymmetric directions, which is useful as a stricter filter.

## Provenance

This setup is adapted from:

- Adaptyv Bio's Nipah workflow: [adaptyvbio/nipah_ipsae_pipeline](https://github.com/adaptyvbio/nipah_ipsae_pipeline)
- The bundled scorer script: [DunbrackLab/IPSAE](https://github.com/DunbrackLab/IPSAE)
- Boltz prediction docs: [jwohlwend/boltz prediction docs](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)
