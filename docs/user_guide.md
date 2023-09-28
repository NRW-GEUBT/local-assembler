# Guide for users

## Installation

### Conda

Install conda from any distribution, i.e. miniconda.
You can follow the setup guide from the [Bioconda team](https://bioconda.github.io/).

We advise installing the mamaba solver in the base environement to speed up
environments creation.

```bash
conda install mamba -n base -c conda-forge
```

Note that you might needto set thechannel priority to 'felxible' to be able to install
the required environements:

```bash
conda config --set channel_priority flexible
```

### Run environment

The running environement simply required a recent python version (>= 3.9) and snakemake.
If you followed the steps above just run:

```bash
mamba create -n snakemake snakemake
```

### Install module and databases

Download the [latest realease](https://github.com/NRW-GEUBT/local-assembler/releases/latest)
and unpack it.

If you're feeling brave, clone the repository form Github:

```bash
git clone https://github.com/NRW-GEUBT/local-assembler
```

Most software and databases dependencies will be installed during the first run.

## Configuration

The configuaration can be defined in two ways:

- either edit and locally save the `config/config.yaml` files and provide its path
  to the snakemake command with the `--configfile` argument

- or provide the parameters directly to the snakemake command with
  `--config <ARGNAME>=<VALUE>`

### User defined parameters

Following arguments must be provided for each run:

| Parameter | Type | Description |
| --- | --- | --- |
| `workdir` | path-like string | Path to the ouptut directory |
| `fastq_folder` | path-like string | Path to the folder containing thesequencing raw-data |
| `metadata` | path-like string | Path to the metadata table intab-delimited format |

### Optional parameters

Following parameters are optional and will revert to defaults if not set:

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `max_threads_per_job` | integer | 1 | Max number of threads assigned to a single job |
| `fastq_naming` | string | `illumina` | FASTQ naming schema(choose from 'illumina', 'ncbi', or 'flex') |

## Usage

The workflow can be started with:

```bash
snakemake --use-conda --conda-prefix <PATH TO CONDA ENVS> --configfile <PATH TO CONFIG> --cores <N>
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more information.

## Input formats

### Sequence files

Raw sequencing data should be provided as two fastq.gz files for each samples (paired-end illumina sequencing only,not interleaved).
File naming should explicitely follow one the 4 schemes:
- illumina: `{samplename}_S*_ L*_R{1,2}_001.fastq.gz`. This is default illumina naming scheme
- ncbi: `{samplename}_{1,2}.fastq.gz`
- flex: `{samplename}*_R{1,2}*.fastq.gz`. The sample name is cut after the first \"_\".

This naming convention is used to generate a sample sheet with the AQUAMIS `create_sampleSheet` helper.
See the documentation directly in the AQUAMIS repository.

## Metadata

Metadata hould be provided as a table ina a flat text format, using tab separators.
The table must contain following fields, not all of which must be filled:

| Field  | Description | Required |
|---|---|---|
| isolate_id | Unique data identification in database | required |
| sample_id | Unique sample identification for epidemiological analysis | required |
| organims | Sample categorization for downstream analysis | required |
| isolate_name_alt | Laboratory internal name | optional |
| isolation_org | Contact organisation for isolate | required |
| sequencing_org | Contact organisation for sequencing | required |
| extraction_method | Extraction method name | optional |
| library_method | Library preparation name | optional |
| sequencing_instrument | Sequencing instrument name | optional |
| bioinformatics_org | Contact organisation fro bionformatic analysis | required |
| third_party_flag | Does this isolate belong to a thrid party organisation? | required |
| thirs_party_owner | Name of the isolate own, required if third_party_flag is TRUE | required |
| sample_type | Sample category | required |
| collection_date | Date ofsample collection | required |
| collection_municipality | Place of sample collection | required |
| collection_country | Country of sample origin | required |
| collection_cause | Reason for sample collection | required |
| collected_by | Local contact for collection | required |
| manufacturer | Production or husbandery firm or facility | optional |
| designation | Product category | required |
| manufacturer_type | Type of production or husbandery facility | required |
| sample_description | Description of the sample | required |
| lot_number | Lot / charge number or stable deignation, etc… | optional |

## Results

### Report

The HTML AQUAMIS report is located under `aquamis/reports/assembly_report.html`.
It contains all QC informations for each sample.

### Assemblies

All assemblies are located under `aquamis/Assembly/assembly`.
This includes assemblies that do not satisfy QC.

### Geübt-ready data

The `geuebt_export` folder contains the data to be sent to the Geübt-core workflow.
This includes all assemblies passing QC and a table containing formatted metadata.
