# BALF EV proteomics: Figure 3 revision

This self-contained analysis folder addresses Reviewer #3, Comment #4. It regenerates panels A–D, with all 1,116 protein-level observations tested in each contrast, from recovered normalized BALF data. It does not change pathway panels E–G or reconstruct upstream PSM processing. The source notebooks in the parent repository should be retained unmodified.

## Run locally

Open a terminal in this folder (the folder containing the notebook).

```bash
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
jupyter notebook BALF_Figure3_revision.ipynb
```

Run all cells. Alternatively, run `python make_volcano.py`. Tested with Python 3.13 and the pinned scientific Python packages in `requirements.txt`. No network access is needed after installation. Generated outputs are written into `figures/` and `results/`. GitHub can preview the executed notebook but does not run it.

## Contents and provenance

- `BALF_Figure3_revision.ipynb`: annotated analysis, validation, plot generation, and preview.
- `make_volcano.py`: equivalent non-interactive analysis script (run from this directory).
- `data/BALF_Normalized_reconstructed.csv`: normalized BALF data with original mouse identifiers; reconstructed from the BALF PSM export using the historical workflow.
- `data/BALF_Normalized_named_reconstructed.csv`: same values, labeled with original notebook group columns.
- `data/sample_map.csv`: mapping of mouse identifiers to group matrix columns, based on the matching column order and the author-confirmed treatment mapping.
- `data/Annotation (2).csv`: original 704 complete-case protein records, used for numerical validation.
- `data/grp* pi` comparison CSVs (`grp21pi (2).csv`, etc.): original π-selected exports, used to verify the selected protein sets and archived p-values, log2 fold changes, and π scores, and to retrieve original labels.
- `figures/`: combined figure and individual panels in vector PDF/SVG and 600-dpi PNG, plus a smaller preview.
- `results/`: all 1,116 proteins per comparison, including group sample counts, means, fold changes, p-values, π scores, selection flags, reference-only BH/Bonferroni values, and label flags. `comparison_summary.csv` summarizes counts. These are working BALF tables, not the final combined serum/BALF supplement.

## Analysis and design decisions

Welch two-sided t-tests use available non-imputed normalized abundance measurements; fold changes use ratios of group arithmetic means. Positive log2 fold change means higher abundance in the first-named group. Colors use |log2 FC| >= 0.5 and nominal p < 0.05. Historical π selection uses |π| > 1.0 and is verified separately. Reference adjusted p-values are not used for figure colors or historical selections. These are exploratory thresholds, not claims of multiple-testing-adjusted significance. A maximum of four labels per direction are selected by smallest nominal p-value among colored points. Gene names are incomplete (the full original annotation table was not recovered); use accession IDs where original labels are unavailable. Shared axes and the lighter blue/orange palette are new presentation choices.

## Repository integration

Place the entire `BALF_Figure3_revision/` folder inside the existing HarryUCLA repository. Keep the project private or confirm coauthor/journal approval before making unpublished data or figures public. Do not overwrite the original notebooks. The raw PSM export and upstream normalization notebook are not duplicated in this package.
