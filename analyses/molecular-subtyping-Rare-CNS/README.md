# Molecular Subtyping of Rare CNS

**Module author:** Aylar Babaei, Jo Lynne Rokita

## Overview

This directory contains scripts and outputs for the Rare CNS molecular subtyping workflow

The goal of this analysis is to identify Rare CNS tumors for downstream molecular subtyping. 

Related issue: [#114](https://github.com/rokitalab/OpenPedCan-Project-CNH/issues/114)

## Workflow steps

The Rare CNS subtyping workflow is organized into the following scripts:

- `00-rare-cns-select-tumors.R`

## Current PR scope

This PR includes only:

- `00-rare-cns-select-tumors.R`
- `rare-cns-subset/rare_cns_subtyping_path_dx_strings.json`
- this README

Downstream steps will be added in stacked follow-up pull requests.

## Step 00

### Purpose

`00-rare-cns-select-tumors.R` identifies candidate Rare CNS samples based on methylation subtype criteria and generates the JSON file for selection in downstream subtyping steps.

### Outputs

- `rare-cns-subset/rare_cns_subtyping_path_dx_strings.json`

### Usage

From within this directory:

```sh
Rscript --vanilla 00-rare-cns-select-tumors.R