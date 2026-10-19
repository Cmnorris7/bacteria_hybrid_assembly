# Public Repository Review and TODO

## Hiring-manager verdict

As currently presented on GitHub, this is not yet a strong public portfolio repository.

The code on `nextflow-refactor` is promising and demonstrates useful production-bioinformatics skills. Unfortunately, GitHub's default branch is still `master`, so a reviewer sees the older shell pipeline rather than the substantial Nextflow work.

Current assessment:

- Public presentation today: 3/10
- `nextflow-refactor` implementation: 6.5/10
- Potential after cleanup and validation: 8/10

## What gives me confidence

- [x] Nextflow DSL2 with sensible workflow, subworkflow, and module separation.
- [x] Local and SLURM execution profiles.
- [x] Versioned containers, with custom module images pinned by digest.
- [x] Multi-sample CSV input rather than filename discovery.
- [x] Reuse of nf-core modules instead of unnecessarily reinventing wrappers.
- [x] Parameter schema, resource labels, reports, trace, timeline, and DAG output.
- [x] A containerized stub test and GitHub Actions workflow.
- [x] A passing remote CI run: https://github.com/Cmnorris7/bacteria_hybrid_assembly/actions/runs/39186877612
- [x] README documentation for inputs, outputs, execution, profiles, and architecture.
- [x] No hard-coded credentials found during review.

These are valuable signals for production bioinformatics and workflow-engineering roles.

## Critical issues

### 1. Put the Nextflow implementation on the default branch

GitHub's default branch is still `master`. A hiring manager opening the repository currently sees:

- An imperative shell pipeline.
- No CI, parameter schema, or tests.
- Destructive input handling in `run_bact_assembly.sh`, which moves FASTQ files.
- Unquoted shell variables and fragile filename discovery using `ls`.
- An installer that clones the current Autocycler HEAD rather than a pinned release.
- A README ending with "Future Work — Add test files."

This hides nearly all the evidence that would help a candidacy.

- [ ] Merge `nextflow-refactor` into the default branch.
- [ ] Make the cleaned Nextflow branch the GitHub default branch.
- [ ] Confirm that the CI workflow passes on the new default branch.
- [ ] Update any repository links or badges that still reference `master`.

### 2. Add scientific validation, not only workflow stubs

`tests/smoke.sh` uses Nextflow stub mode. This demonstrates that mocked process invocations connect, but it does not establish that:

- A real bacterial isolate completes successfully.
- The resulting assembly is accurate.
- Polishing improves the assembly.
- Plasmids are retained.
- Multiple samples remain isolated correctly.
- Local and SLURM profiles produce comparable results.

- [ ] Select a small public isolate with a trusted reference assembly.
- [ ] Document its accession, read chemistry, coverage, and expected genome size.
- [ ] Add a reproducible real-data integration test or validation command.
- [ ] Record expected output checksums or bounded acceptance criteria.
- [ ] Demonstrate multi-sample execution with at least two isolates.
- [ ] Document whether local and SLURM runs produce equivalent biological outputs.

### 3. Strengthen final assembly QC

The current final assessment primarily maps Nanopore reads back to the assembly and runs a custom shell script. It lacks standard evidence such as:

- QUAST metrics and reference comparison.
- BUSCO completeness.
- CheckM2 completeness and contamination.
- Short-read-based consensus quality or residual variant counts.
- Expected chromosome and plasmid recovery.
- Comparison against a known truth assembly.

The custom circularity score is especially risky. Its thresholds are hand-built and uncited, and it could communicate more confidence than the evidence supports.

- [ ] Add QUAST to the final QC workflow.
- [ ] Consider BUSCO and/or CheckM2, with documented database requirements.
- [ ] Add a consensus-quality or residual-error assessment using short reads.
- [ ] Report chromosome and plasmid recovery for the validation dataset.
- [ ] Cite or validate every custom QC threshold.
- [ ] Label the custom circularity analysis as exploratory until it is validated.
- [ ] Prefer established circularity evidence or tooling where possible.
- [ ] Add tests for `assemblyQC_coverage_stats.sh` using controlled BAM fixtures.

### 4. Enforce input and parameter validation

The README says sample IDs must be unique, but `input_qc.nf` does not enforce uniqueness. The parameter schema also does not appear to be connected to runtime validation.

Current gaps include:

- `help` is defined but not implemented.
- `additionalProperties` is `true`.
- The three-distinct-assemblers requirement is not encoded in the schema.
- CSV structure, duplicate samples, file extensions, and mate consistency are not comprehensively validated.

- [ ] Add formal samplesheet validation.
- [ ] Reject duplicate sample IDs with a clear error.
- [ ] Validate required columns before constructing channels.
- [ ] Validate that all read paths exist and have supported compression/extensions.
- [ ] Validate paired-read consistency and prevent R1/R2 reuse across samples.
- [ ] Connect `nextflow_schema.json` to runtime validation or document its exact purpose.
- [ ] Implement `--help` or remove it from the schema.
- [ ] Reject unknown parameters where practical.
- [ ] Enforce at least three distinct supported assemblers.
- [ ] Add negative tests for malformed samplesheets and invalid parameters.

### 5. Finish the reproducibility story

Containers are a strong start, but several moving parts remain mutable or unpublished.

- [ ] Pin the Nextflow version used by CI.
- [ ] Pin GitHub Actions to immutable commit SHAs.
- [ ] Aggregate module version topics into a published `software_versions.yml` or equivalent.
- [ ] Include the Nextflow version, pipeline commit, command, and profile in the results.
- [ ] Create a tagged release for the first validated Nextflow version.
- [ ] Add a changelog or release notes.
- [ ] Document how vendored nf-core modules are updated and reviewed.

### 6. Remove the competing legacy interface

The README says `main.nf` is supported, but the repository root still contains `run_bact_assembly.sh` and `setup_environments.sh`, while `src/` contains the old implementation. This makes it unclear which interface users should trust.

- [ ] Move the legacy shell pipeline to an archival branch or release tag.
- [ ] Remove unsupported root-level entry points from the primary branch.
- [ ] Keep only scripts still required by the Nextflow workflow.
- [ ] Clearly mark any retained legacy code as unsupported and non-production.
- [ ] Ensure no supported command moves or modifies raw input files.

### 7. Improve repository and commit hygiene

The Nextflow conversion is currently one very large commit:

- 119 files changed.
- 10,445 additions.
- 81 deletions.

Earlier history also contains messages such as "need to fix," "need to test," and "WIP." This is normal during development, but it is not ideal as the visible history of a showcase repository.

- [ ] Use focused, reviewable commits for future work.
- [ ] Use descriptive commit messages that explain intent and verification.
- [ ] Merge the refactor through a documented pull request if practical.
- [ ] Add release notes summarizing the architecture change.

### 8. Add public-project metadata

The GitHub repository currently has no description, topics, license, release, or citation information.

- [ ] Add a clear GitHub repository description.
- [ ] Add topics such as `bioinformatics`, `nextflow`, `bacterial-genomics`, `hybrid-assembly`, `nanopore`, and `illumina`.
- [ ] Choose and add a license.
- [ ] Add `CITATION.cff` if the workflow should be cited.
- [ ] Add a CI status badge to the README.
- [ ] Add a version/release badge after publishing the first release.
- [ ] Add one representative DAG, report screenshot, or results figure.
- [ ] Add an example results summary from the validation dataset.

## Recommended order of work

### Priority 1: Stop underselling the project

- [ ] Make the Nextflow implementation the default branch.
- [ ] Remove or archive the legacy shell implementation.
- [ ] Add a license, repository description, and topics.
- [ ] Add a CI badge and a prominent project-status statement.

### Priority 2: Demonstrate biological trustworthiness

- [ ] Add one reproducible real-data validation case.
- [ ] Publish expected QUAST, BUSCO, and/or CheckM2 results.
- [ ] Demonstrate that polishing improves a measurable metric.
- [ ] State limitations clearly, especially around circularity and plasmid recovery.

### Priority 3: Demonstrate production discipline

- [ ] Enforce samplesheet and parameter validation.
- [ ] Aggregate and publish tool versions.
- [ ] Pin Nextflow and CI actions.
- [ ] Add nf-test or focused module/subworkflow tests.
- [ ] Create a tagged, documented release.

## Hiring perspective

For a bioinformatics workflow engineer position, the Nextflow branch is strong enough to prompt a discussion about design choices. It demonstrates real domain understanding and useful infrastructure instincts.

For a senior production or Staff/Principal-level position, the repository alone is not yet convincing. The missing scientific validation, weak release discipline, and default-branch presentation make it look like a capable internal pipeline that has not completed the transition into a trustworthy public software product.

The underlying work is worth developing. The current public presentation materially undersells it.
