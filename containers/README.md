# Reusable GeneInIt container builds

These recipes build the two GeneInIt runtime layers without project-specific
paths:

1. `predictor/` builds GeneMark-ETP 1.02 and AUGUSTUS 3.5.0.
2. `allinone/` adds GeneInIt, Snakemake, RNA evidence tools, and InGenAnnot.

The build scripts derive the GeneInIt source root from their own locations.
Output, licensed inputs, and image names are supplied through environment
variables. Neither Docker image nor SIF contains a GeneMark license key.

## Prerequisites

- Docker with BuildKit
- Apptainer
- A locally obtained `GeneMark-ETP.tar.gz`
- A valid GeneMark key
- The GeneMark standalone archive containing its bundled smoke-test fixture

GeneMark software and license files remain subject to their upstream license.
Do not commit or redistribute them with GeneInIt.

## Build the predictor layer

```bash
GENEMARK_ETP_ARCHIVE=/path/to/GeneMark-ETP.tar.gz \
GENEMARK_KEY_PATH=/path/to/gm_key \
GENEMARK_STANDALONE_ARCHIVE=/path/to/gmes_linux_64_4.tar.gz \
GENEINIT_OUTPUT_DIR=/path/to/output \
containers/predictor/build_sif.sh
```

Optional predictor settings are `GENEINIT_PREDICTOR_IMAGE`,
`GENEINIT_SIF_NAME`, `GENEINIT_SOURCE_CACHE_DIR`, `AUGUSTUS_ARCHIVE`, and
`AUGUSTUS_URL`.

## Build the all-in-one layer

Build or otherwise load the predictor Docker image first, then run:

```bash
GENEMARK_KEY_PATH=/path/to/gm_key \
GENEMARK_STANDALONE_ARCHIVE=/path/to/gmes_linux_64_4.tar.gz \
GENEINIT_OUTPUT_DIR=/path/to/output \
containers/allinone/build_sif.sh
```

Optional all-in-one settings are `GENEINIT_PREDICTOR_IMAGE`,
`GENEINIT_ALLINONE_IMAGE`, `GENEINIT_BUILD_REVISION`, and
`GENEINIT_SIF_NAME`.

Both scripts fail if the requested final SIF already exists. After conversion,
they verify pinned versions, all GeneMark-ETP Perl entry points, required CPAN
modules, and a licensed GeneMark core prediction. The key is copied with mode
`0600` into a temporary Apptainer home only and is removed when the build exits.

The all-in-one image deliberately uses separate alignment and StringTie Conda
environments. This preserves samtools 1.24 while avoiding its incompatible
HTSlib constraint with StringTie 3.0.3.

## Run with an explicit CPU budget

```bash
apptainer exec --cleanenv --bind /data:/data /path/to/geneinit.sif \
  geneinit-run --threads 34 --configfile /data/config.yaml
```

`--threads` is the total budget for the workflow, including concurrent jobs.
The older numeric `--cores` / `-c` options remain accepted. Execution without a
budget is rejected; `--help` and `--version` remain available without one.
The launcher forwards the budget to Snakemake and restricts inherited Linux CPU
affinity. The image includes sequence-parallel AUGUSTUS prediction and its
official output merger. See the main README for per-stage behavior.
