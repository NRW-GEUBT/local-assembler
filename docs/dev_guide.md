# Guide for developers

## Changing quality thresholds

To change quality thresholds, modify the JSON file located in
`workflow/schema/AQUAMIS_threesholds.json`.

The file provided in local-assembler is a scaled-down, modified version of the
[AQUAMIS definition](https://gitlab.com/bfr_bioinformatics/AQUAMIS/-/blob/master/resources/AQUAMIS_thresholds.json?ref_type=heads).

The modifications are:
- removed unused organisms
- modified some of the thresholds to follow those defined in official normative methods.

## Changing AQUAMIS configuration

The parameters for AQUAMIS execution are currently fixed, in order to limit user burden as well as to
ensure consistent execution.
These parameters can bechanged directly in the snakemake rule executing AQUAMIS in the following file:
`workflow/rules/aquamis.smk`.
