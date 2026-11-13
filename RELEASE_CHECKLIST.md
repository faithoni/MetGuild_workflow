# Release acceptance checklist

Candidate: **0.1.0-rc.1**. Canonical version: [VERSION](VERSION).
This checklist records remaining acceptance work; unchecked items are not claimed as passed. Benchmarking is outside this release-preparation scope.

## Fresh environment and Bash launcher

Run from a checkout of the exact candidate commit in a Bash-capable environment with Conda initialized. Record the OS, Bash/Conda versions, commit SHA, resolved package versions, commands, logs, and exit codes. Use a new environment name and fresh output paths if these already exist.

```bash
git rev-parse HEAD
cat VERSION
bash --version
conda --version
conda env create -n metguild_rc_010rc1 -f environment.yml
conda activate metguild_rc_010rc1
conda list
bash -n run_metguild_workflow.sh
bash run_metguild_workflow.sh --help

bash run_metguild_workflow.sh Example_dataset/binary_matrix.csv 15 \
  --matrix-type binary --taxonomy Example_taxonomy.csv \
  --output-dir rc_binary_with_taxonomy
echo "Binary exit code: $?"

bash run_metguild_workflow.sh Example_dataset/proportional_matrix.tsv 15 \
  --matrix-type proportion --taxonomy Example_taxonomy.csv \
  --output-dir rc_proportion_with_taxonomy
echo "Proportional exit code: $?"

bash run_metguild_workflow.sh Example_dataset/binary_matrix.csv 15 \
  --matrix-type binary --output-dir rc_binary_without_taxonomy
echo "Without-taxonomy exit code: $?"
```

- [ ] Fresh environment creation, Bash syntax check, and help command succeed.
- [ ] All three launcher runs return exit code 0 and reach Step 9.
- [ ] Both clustering methods produce outputs consistent with the candidate-k manifest.
- [ ] Network edge endpoints match node tables; manifest row counts and checksums agree with files.
- [ ] Taxonomy appears in cluster attributes when supplied; the run without taxonomy succeeds.
- [ ] SQLite, SQL, and HTML exports exist in each run's 08_workflow_provenance directory.
- [ ] SQLite integrity/foreign-key checks pass, and the SQL dump recreates the database (allow normal floating-point text round-trip differences).
- [ ] Existing output directories are rejected without overwriting results; an invalid matrix type and conflicting stop flags return nonzero.
- [ ] Stop modes finish at their documented stage without claiming full-run provenance export.

## Real-browser acceptance

Open both taxonomy-enabled reports directly from their complete run directories:

- rc_binary_with_taxonomy/08_workflow_provenance/binary_matrix_workflow_provenance.html
- rc_proportion_with_taxonomy/08_workflow_provenance/proportional_matrix_workflow_provenance.html

Record browser name/version and any console errors.

- [ ] Reports render and graph zoom works.
- [ ] Search and step, result-type, method, and candidate-k filters work independently and in combination.
- [ ] Method/k filters retain shared upstream sources and exclude unrelated solutions.
- [ ] Selecting a node shows its correct path, metadata, preview, and upstream/downstream relationships.
- [ ] Each filter refreshes details for a still-visible selection or clears details for a hidden selection.
- [ ] The hierarchical candidate-k figure connects to hierarchical solutions.
- [ ] Table and figure previews display, and representative CSV/TSV, figure, and interactive-plot links open the correct existing files.
- [ ] Clearing filters restores the map without stale details.

## Promotion record

- [ ] Record the tested commit SHA and attach or link acceptance evidence.
- [ ] Resolve any blocking failures and repeat affected checks.
- [ ] Confirm VERSION, README, and Releases.md agree on the intended release version.
- [ ] When publishing this candidate, use tag v0.1.0-rc.1 and mark the GitHub release as a prerelease.
- [ ] Before a stable release, complete outstanding acceptance checks and update version and release notes.

Tag creation and publication are separate from preparing these files.
