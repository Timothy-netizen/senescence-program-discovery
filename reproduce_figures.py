"""Regenerate saved-result figures without downloading data or fitting models."""

import argparse
import ast
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
NOTEBOOKS = tuple(
    f"{name}_results/{name}_notebook.ipynb"
    for name in ("geometry", "bars", "evidence", "biological")
)


def run_notebook(relative_path):
    from IPython.display import display

    notebook = json.loads((ROOT / relative_path).read_text(encoding="utf-8"))
    sources = ["".join(cell["source"]) for cell in notebook["cells"]
               if cell["cell_type"] == "code"]
    assignments = [node for source in sources for node in ast.walk(ast.parse(source))
                   if isinstance(node, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == "RUN_FITS"
                           for t in node.targets)]
    if len(assignments) != 1 or ast.literal_eval(assignments[0].value) is not False:
        raise RuntimeError(f"Set RUN_FITS = False in {relative_path} before regenerating figures.")
    namespace = {"__name__": "__main__", "display": display}
    for index, source in enumerate(sources):
        exec(compile(source, f"{relative_path}:cell{index + 1}", "exec"), namespace)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--notebook", choices=NOTEBOOKS, help=argparse.SUPPRESS)
    args = parser.parse_args()
    os.chdir(ROOT)
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    if args.notebook:
        run_notebook(args.notebook)
        return
    for notebook in NOTEBOOKS:
        print(f"Regenerating {notebook}", flush=True)
        subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()),
                        "--notebook", notebook], check=True)
    subprocess.run([sys.executable, "-B", "senescence_results/make_thesis_figures.py"], check=True)
    print("Done. Figures are in the experiment folders and senescence_results/thesis_figures/.")


if __name__ == "__main__":
    main()
