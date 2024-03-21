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

Note that the workflow will expect fastq to be named after either the value of the `isolate_name_alt` if it is filled
or the value of the `isolate_id` field if `isolate_name_alt` is empty. If a fastq has a name that is not in these fields,
the workflow will stop with an error.

## Metadata

Metadata should be provided as a table in a flat text format, using tab separators.
The table must contain following fields, not all of which must be filled:

| Field name             | Definition                                                                                            | Required? | Allowed Values                                                                                 |
| ---------------------- | ----------------------------------------------------------------------------------------------------- | --------- | ---------------------------------------------------------------------------------------------- |
| isolate_id             | Eindeutiger Bezeichner der Isolat                                                                     | required  | free text                                                                                      |
| sample_id              | Eindeutiger Bezeichner der Probe                                                                      | required  | free text                                                                                      |
| alt_isolate_id         | Alternativer Bezeichnung der Siolat (z.B. Laborrinternenummer)                                        | optional  | free text                                                                                      |
| organism               | Erwarteter Organismus                                                                                 | required  | oneof "Listeria monocytogenes","Salmonella enterica", "Escherichia coli", "Campylobacter spp." |
| isolation_org          | Identifizierungsname des für die Isolation verantwortlichen Labor                                     | required  | oneof "OWL", "RRW", "WFL", "RLD","MEL", "other"                                                |
| sequencing_org         | Identifizierungsname des für die Sequenzierung verantwortlichen Labor                                 | required  | oneof "OWL", "RRW", "WFL", "RLD","MEL", "other"                                                |
| bioinformatics_org     | Identifizierungsname des für die Primär Analyse verantwortlichen Labor                                | required  | oneof "OWL", "RRW", "WFL", "RLD","MEL", "other"                                                |
| third_party_owner      | Ggfs. Name der Dritteigentümer (z.B. bei Humanem Isolaten oder RV)                                    | optional  | free text                                                                                      |
| extraction_method      | Name des Kits oder Methode für die DNA Extraktion                                                     | optional  | free text                                                                                      |
| library_kit            | Name des Kits für die DNA-Library Erstellung                                                          | optional  | free text                                                                                      |
| sequencing_kit         | Name des Kits für die DNA Sequenzierung                                                               | optional  | free text                                                                                      |
| sequencing_instrument  | Model des Gerät für die DNA Sequenzierung                                                             | optional  | free text                                                                                      |
| assembly_method        | Name und Version des verwendeten Software oder Pipeline für die Assembly                              | optional  | free text                                                                                      |
| sample_type            | Art der Probe (z.B. Lebensmittel, Umfeldprobe, veterinärmedizinisches Material, usw.)                 | required  | oneof "Lebensmittel", "Futtermittel", "Umfeld", "Tiergesundheit", "Human", "unknown"           |
| collection_date        | Datum der Probenahme                                                                                  | required  | ISOdate                                                                                        |
| customer               | Einsendende KOB                                                                                       | required  | Kodierung (2-5 buchstaben)                                                                     |
| manufacturer           | Hersteller der Probe                                                                                  | required  | Freitext OR  "unknown"                                                                         |
| collection_place       | Entnahmeort der Probe                                                                                 | required  | AVV / ADV Text                                                                                 |
| collection_place_code  | Entnahmeort der Probe                                                                                 | required  | AVV / ADV Code                                                                                 |
| description            | Bezeichnung laut Probeentnahmeschein bzw. bei Tupfern Bezeichnung der Entnahmeortes und ggfs. Tierart | required  | free text                                                                                      |
| manufacturer_type      | Einordnung des Entnahmeortes nach AVV-Katalog                                                         | optional  | AVV Text                                                                                       |
| manufacturer_type_code | Einordnung des Entnahmeortes nach AVV-Katalog                                                         | optional  | AVV Code                                                                                       |
| matrix                 | Einordnung der Probenmatrix nach AVV-Katalog                                                          | optional  | AVV/TSN Text                                                                                   |
| matrix_code            | Einordnung der Probenmatrix nach AVV-Katalog                                                          | optional  | AVV/TSN Code                                                                                   |
| collection_cause       | Grund für die Probenahme nach AVV-Katalog                                                             | optional  | AVV Text oder Freitext                                                                         |
| collection_cause_code  | Grund für die Probenahme nach AVV-Katalog                                                             | optional  | AVV Code                                                                                       |
| lot_number             | Los-Nr./Chargen-Nr. des Herstellers                                                                   | optional  | free text                                                                                      |

Note that in `local-assembler` no in-depth check of the provided metadata is performed.
This will however happen during the submission to geuebt-core.

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
