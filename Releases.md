# Release notes

## 0.1.0 — stable release

Prepared on 2026-11-19. Promotes the workflow source from 0.1.0-rc.2 to 0.1.0. Analysis algorithms, statistical defaults, and output formats are unchanged.

### Changed

Relative input, taxonomy, and output paths now resolve from the directory where the command is invoked. Automatically numbered output directories are also created there. Bundled scripts continue to resolve relative to the launcher. Running the documented examples from the repository directory retains the same behavior.

### Included functionality

- Complete workflow for binary presence/absence and proportional genome-by-feature matrices, with optional taxonomy annotation.
- Candidate-k selection, classical-MDS/KMeans and average-linkage hierarchical clustering, statistical screening, and filtered Spearman association networks.
- Network edge lists and node attributes with solution-specific cluster identifiers.
- SQLite provenance database, SQL reconstruction dump, and interactive HTML workflow map with search, filters, previews, and result links.
- Bundled binary and proportional example datasets with matching taxonomy.

### Fixes carried forward

- pandas 3 compatibility for correlation-matrix updates.
- Preservation of string genome identifiers and numeric-looking feature names.
- Distinct feature and cluster-node identifiers in network exports.
- Hierarchical candidate-k figure connections, details-panel refresh and clearing when filters change, and corrected HTML text encoding.

### Usage

Follow the [README quick start](README.md#quick-start). Keep the complete output directory when sharing HTML reports so previews and result links remain accessible.

This version change does not assert additional test coverage; the [release checklist](RELEASE_CHECKLIST.md) retains the recorded status. Tag creation and GitHub publication are separate from this source preparation.

## 0.1.0-rc.2 — release candidate

Prepared on 2026-11-16.

### Fixed

- Restored explicit lineage from the cleaned matrix to the hierarchical candidate-k figure, and from that figure to hierarchical clustering solutions.
- Restored details-panel updates when search or filters change: visible selections refresh, and hidden selections clear.
- Corrected corrupted ellipsis, bullet, and minus characters in HTML previews, graph labels, and zoom controls.

### Release metadata

Updated VERSION, README, and the release acceptance checklist to 0.1.0-rc.2. This candidate restores the HTML-map behavior documented for rc.1; analysis algorithms and statistical defaults are unchanged. Preparing these files does not publish a Git tag or GitHub release.

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
