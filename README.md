# Senescence program discovery

Code and saved results accompanying Timothy Lauren Nijhuis's MSc thesis on sparse dictionary learning and senescence-associated expression programs. This repository provides a quick route to regenerate figures and the code for repeating the experiments and atlas analysis.

## Contents

| Location | Contents |
| --- | --- |
| `margin_dictionary_learning.py` | Dictionary-learning implementation |
| `geometry_results/`, `bars_results/` | Geometric and bars experiments: notebooks and saved scores |
| `evidence_results/`, `biological_results/` | Repeated recovery experiments: notebooks and saved scores |
| `senescence_results/` | Atlas preprocessing, calibration, fitting and analysis code; calibration and fitting summaries |
| `senescence_results/analysis/` | Program rankings, enrichment results and recorded molecular interpretations |
| `senescence_results/references/` | Frozen biological reference inputs |

Thesis text, bibliography, large expression matrices, fitted atlas files, checkpoints and temporary downloads are omitted. Generated figures are also omitted because the included results reproduce them.

## Setup

Use Python 3.13 in a virtual environment. From the repository root:

```sh
python -m venv .venv
```

Activate it with `.venv\Scripts\Activate.ps1` in Windows PowerShell, or `source .venv/bin/activate` on Linux/macOS, then install the recorded local dependencies:

```sh
python -m pip install -r requirements.txt
```

Open notebooks in a Jupyter-compatible editor using this environment. Full-atlas calibration and fitting instead target Google Colab with an A100 40-GB GPU, CUDA PyTorch and Triton. The local requirements file does not replace that GPU environment; recorded calibration versions are in `senescence_results/calibration_metadata.json`.

## Regenerate figures from saved results

Run:

```sh
python reproduce_figures.py
```

This executes the four experiment notebooks with their saved-results setting, `RUN_FITS = False`, then runs `senescence_results/make_thesis_figures.py`. It requires no atlas download or model fitting. Experiment PDFs appear in their respective results folders; atlas figures appear in `senescence_results/thesis_figures/`. The saved CSV files remain unchanged.

The notebooks also display numerical summaries. The atlas ranking is in `analysis/screen/program_scores.csv`; the seven selected candidates and their molecular evidence are under `analysis/top_7_4fa212641968/`, relative to `senescence_results/`.

## Repeat fitting and analysis

Work in a separate clone when refitting: experiment notebooks replace their score CSVs. Set `RUN_FITS = True` and keep `RESUME_FITS = False` to repeat calibration and evaluation from the beginning. Seeds, grids and numerical settings are defined near the beginning.

For the atlas, follow this order:

1. Run `python senescence_results/prepare_atlas.py`. It downloads the pinned Allen atlas release, checks publisher checksums and creates `senescence_atlas_processed.h5` beside the script. Downloads and preprocessing resume after interruption. Allow space for approximately 14 GB of source downloads, a 14 GB processed file and working files.
2. Upload the processed file to Google Drive and the shared Python module to Colab. Run `senescence_calibration.ipynb`, adjusting its Drive paths and choosing a fresh `OUTPUT_DIR`.
3. Run `senescence_fit.ipynb`. Set `CALIBRATION_DIR` to that new calibration directory, which must contain its CSV, `experiment.npz` and generated checkpoints. Choose a fresh fit output directory. The saved calibration CSV alone cannot initialise this fit.
4. Copy the resulting `atlas_fit.h5` into `senescence_results/` in your working clone. Run both cells of `senescence_analysis.ipynb` in order. It uses the included frozen biological inputs. New fits can change program IDs and require fresh molecular interpretation; the supplied figure script targets the saved thesis results.

## Reference data

The included reference inputs contain Reactome 97 pathways, the mouse SenMayo panel and gene-identifier mappings.
