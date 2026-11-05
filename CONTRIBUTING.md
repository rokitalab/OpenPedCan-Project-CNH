# Contribution guidelines for OpenPedCan-Project

OpenPedCan is described in [Geng et al., _GigaScience_ (2025)](https://doi.org/10.1093/gigascience/giaf093), but the project is a living repository: the data releases, analyses and documentation continue to be updated and extended.
Anyone is welcome to contribute, whether by asking a question, reporting a data problem, proposing an analysis or submitting a pull request.

Thank you for your interest in contributing to OpenPedCan.
Please read the [README](./README.md) first, in particular the sections on the data release process and [How to Add an Analysis](./README.md#how-to-add-an-analysis).
Data files and formats are described in [`doc/data-formats.md`](./doc/data-formats.md) and [`doc/data-files-description.md`](./doc/data-files-description.md).

## Issues

Use the [issues](https://github.com/rokitalab/OpenPedCan-Project-CNH/issues) to engage with the project.
Please choose the [issue template](./.github/ISSUE_TEMPLATE) that fits, for example:

- `propose-an-analysis`: propose a new analysis module.
- `update-an-analysis`: update an existing analysis module.
- `data-question`: ask about or report a problem with OpenPedCan data.
- `data-release`: track the work for a data release.
- `bix-dev-feature-request`, `bix-ops-harmonization-request`, `bixu-support-pedcbio`: development, harmonization and PedcBioPortal requests.
- `manuscript`: manuscript-related requests.

Before starting an analysis, check the existing issues for a similar one and comment on it to say you plan to work on it.
If there is none, open a new issue describing the scientific goals, methods, input data and timeline.
Post progress updates on the issue, and reference it in your pull request (e.g., `Closes #123`).

## Pull requests

Changes are made by pull request against the `dev` branch of [`rokitalab/OpenPedCan-Project-CNH`](https://github.com/rokitalab/OpenPedCan-Project-CNH).
Do not commit directly to `dev`.

1. Fork the repository (or, if you are a project team member, branch from `dev` in the main repository).
2. Create a branch from `dev` named for the change (e.g., `gene-expression-clustering`).
3. Commit and push your changes, then open a pull request against `dev`.
4. Fill out the [pull request template](./.github/PULL_REQUEST_TEMPLATE.md).
5. Delete the branch once the pull request is merged.

To keep your fork's `dev` current, add the main repository as `upstream`, then:

```
git fetch upstream
git checkout dev
git merge upstream/dev
git push
```

Branches must be up to date with `dev` before merging.
Use the `Update branch` button on the pull request or merge `dev` into your branch; ask a maintainer for help with merge conflicts.

### Size of pull requests

Keep pull requests small and focused on one analysis module or one change, so that reviewers can review the code in detail.
A pull request over roughly 400 lines of code should usually be split into several smaller pull requests.
If you are unsure how to split one up, open a draft pull request and ask for feedback.

### Continuous integration

Pull requests trigger GitHub Actions workflows ([`build.yml`](./.github/workflows/build.yml) and [`run_analysis.yml`](./.github/workflows/run_analysis.yml)) that build the project Docker image and run analysis modules on the CI subset files created by the [`create-subset-files`](./analyses/create-subset-files) module.
Analyses must run in the project Docker container (see [Docker Image](./README.md#docker-image)).
If you add a module that should be tested in CI, add it to the workflow and make sure its inputs are included in the subset files.

## Adding or changing an analysis

Follow [How to Add an Analysis](./README.md#how-to-add-an-analysis) for the folder structure, documentation and script numbering.
In short:

- Put the module in its own folder in [`analyses/`](./analyses), with a `README.md` describing the purpose, methods, inputs, outputs and how to run it, and a `run-*.sh` script that runs it from start to finish.
- Read inputs from the data release files in `data/` and write results to the module's `results/` and `plots/` folders.
- Use relative paths, set seeds for anything stochastic, and comment non-obvious steps.
- Record new software dependencies in the Docker image rather than relying on a local install.
- If the module's output is part of the data release, also update [`doc/data-files-description.md`](./doc/data-files-description.md) and the release notes in [`doc/release-notes.md`](./doc/release-notes.md).
- If the module is part of the data release process, add it to the module tables in the README.

## Peer review

Pull requests are reviewed by the project team ([CODEOWNERS](./.github/CODEOWNERS)) and require at least one approving review before merging.
Reviewers check that:

- The code is correct and follows basic style conventions.
- The results can be reproduced in a reasonable amount of time and are as expected.
- The documentation describes the purpose, methods, results, inputs, outputs and how to run the code.
- Non-obvious procedures are commented.

When reviewing, say what you did (read the code, read the documentation, ran the code, interpreted the results).
Explain why a change is needed and how to make it, and comment on the specific changed lines where possible.
Extensive refactoring that is out of scope for the pull request should go in a new issue.

Do not merge an approved pull request right away.
Coordinate with the maintainers first, since merge order affects other open branches.
