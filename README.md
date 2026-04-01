# High-Throughput Generative-to-DFT Screening

This repository collects scripts for running high-throughput materials screening workflows that connect structure generation, fast prescreening, DFT job orchestration, and post-processing on HPC systems.

The codebase started from electride-focused projects, so some filenames and utilities still use electride-specific naming. The repository itself is broader than that original use case and can be adapted to other screening targets by changing the generation inputs, screening criteria, and analysis steps.

## What Is Here

This is a script-first repository rather than a packaged Python library. Most workflows are organized as small Python and shell tools that can be composed into campaign-specific pipelines.

At a high level, the repository covers:

- generative structure proposal workflows
- fast stability prescreening before DFT
- batch submission and monitoring of DFT jobs
- refinement and follow-up calculations
- result analysis, plotting, and report building

## How To Explore The Repo

Start with these top-level locations:

- `generation/`: structure generation workflows and composition-space exploration scripts
- `prescreen/`: original prescreening workflow
- `prescreen_new/`: newer prescreening variant and filtering utilities
- `origin_flow/`: initial end-to-end DFT workflow scripts, job management, and analysis tools
- `refined_flow/`: refinement-stage workflows and follow-up calculations
- `results_plot/`: plotting and reporting scripts for completed campaigns
- `test/`: small utility and test scripts

Useful top-level files:

- `WORKFLOW_DIAGRAM.md`: high-level workflow sketch
- `create_environment.md`: environment setup notes
- `LICENSE.md`: repository license

## Typical Workflow Shape

Most projects built from this repository follow a pattern like:

1. generate or collect candidate structures
2. prescreen candidates with fast surrogate or filtering steps
3. launch and monitor DFT calculations
4. refine promising candidates
5. analyze, visualize, and export results

You do not need every directory for every campaign. In many cases, a project will use only one generation path, one prescreening path, and one analysis path.

## Practical Notes

- Many scripts are written for HPC and SLURM-based execution.
- The repository contains historical and workflow-specific filenames; treat the directory structure and script behavior as the main guide, not only the names.
- If you are new to the repo, begin with `generation/`, `prescreen_new/`, `origin_flow/`, and `refined_flow/` in that order.

## Current Status

This repository is being reorganized into a more general-purpose screening toolkit. Some scripts remain specialized to earlier electride campaigns, but the overall structure is intended to support broader high-throughput generative-to-DFT studies.
