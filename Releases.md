# Release notes

## 0.1.0-rc.1 — release candidate

Prepared on 2026-11-13. This is a prerelease candidate for validation, not a stable release. The source version is recorded in [VERSION](VERSION); preparation does not create a Git tag or publish a GitHub release.

### Included functionality

- Complete Steps 1–9 workflow for binary and proportional genome-by-feature matrices, with optional taxonomy annotation.
- Candidate-k selection, classical-MDS/KMeans and average-linkage hierarchical clustering, statistical screening, and filtered Spearman association networks.
- Network edge lists and node attributes, with solution-specific cluster identifiers and optional genome–phylum annotations.
- Populated SQLite provenance database, SQL reconstruction dump, and interactive HTML workflow map with filters, details, previews, and relative result links.
- Bundled binary example (127 genomes × 227 features), proportional example (127 × 29 features), and matching taxonomy table.

### Fixes included

- pandas 3 compatibility for correlation-matrix diagonal updates.
- Preservation of string genome identifiers and numeric-looking feature names in downstream outputs.
- Separation of feature names such as C1 from solution-specific cluster nodes during network construction.
- Curated HTML lineage and method/k filtering, including the hierarchical candidate-k plot.
- Details-panel refresh or clearing when filters change the visible selection.

### Validation status

Prior review exercised the Python analysis stages and the runner's embedded network validation and provenance export for both bundled examples with taxonomy and MAX_K=15. Network consistency, SQLite integrity, SQL reconstruction, generated JavaScript syntax, and programmatic map filtering/details behavior were checked.

The reviewed environment was Windows with Python 3.12, NumPy 2.4.2, pandas 3.0.0, SciPy 1.17.0, scikit-learn 1.8.0, Matplotlib 3.10.8, and Plotly 6.5.2. This was an existing Python environment, not a fresh installation from environment.yml. These checks are prior review evidence, not an automated test suite shipped with this candidate.

Full execution through the Bash launcher in a fresh Conda environment and interactive testing of the generated maps in a real browser remain outstanding. Follow [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) and record the candidate commit, environment, logs, and results before stable release promotion. Other allowed Python/dependency versions have not all been tested.

### Usage and limitations

Use the [README quick start](README.md#quick-start) after activating the Conda environment. MAX_K is an exclusive upper bound. Relative input, taxonomy, and output paths resolve against the launcher directory; explicit output directories must be fresh.

A full run requires at least one retained network edge across clustering solutions before provenance export. Keep the complete run directory when opening or sharing an HTML map because figure previews and result links use relative paths.

Results are exploratory associations and require review of statistical diagnostics and domain interpretation. Final manuscript citation, DOI, and remaining manuscript dataset links are pending as described in the README.
