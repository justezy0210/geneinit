from snakemake.utils import validate


validate(config, "config/config.schema.yaml")

include: "workflow/rules/common.smk"


rule all:
    input:
        RELEASE_COMPLETE


include: "workflow/rules/preflight.smk"
include: "workflow/rules/rna.smk"
include: "workflow/rules/prediction.smk"
include: "workflow/rules/evidence.smk"
include: "workflow/rules/models.smk"
include: "workflow/rules/selection.smk"
include: "workflow/rules/coding_models.smk"
include: "workflow/rules/postqc.smk"
include: "workflow/rules/utr_isoforms.smk"
include: "workflow/rules/release.smk"
