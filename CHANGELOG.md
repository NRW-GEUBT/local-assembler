### 1.2.1

- Change upper warning for Listeria SNV count from 3 to 11. This is to reflect experience of the NRL with many cases of Listeria showing high SNV counts.
- Change upper fail for Salmonella assembly length to 6.5 Mb. This reflects experience of the NRL with large Salmonella assemblies.

### 1.2.0

- Add `--fix_fails`to AQUAMIS run options
- Add the possibility of specifying a run name to AQUAMIS with the `run_name` parameter

### 1.1.2

- More ANSI Fix

### 1.1.1

- Fix parsing failure if metadata.tsv was written under windows with ANSI encoding

### 1.1.0

- Calculate and add FASTQ md5 to metadata JSON for Geübt export

### 1.0.4

- Multiple Bugfix NRL export

### 1.0.3

- Fix filtering of samples to export for NRLs

### 1.0.2

- Symlinks are now absolute paths

### 1.0.1

- remove metadata validation, this will take place exclusively at the API level

### 1.0.0

- Update to Geübt Metadata V4
- Export to NRL uses new ADV/AVV template
- Export reads and assemblies are now symlinks

### 0.3.3

Fix NRL xport convenience tool to correclty expor data by organisms

### 0.3.2

Fix consolidation with redundant sample IDs

### 0.3.1

Fix field seq_depth

### 0.3.0

Update Metadtaa to GeÜBT v2
Relaxed validation to accomodate with untypical Sample and isolate naming

### 0.2.4

Modified the naming of FASTQ files for NRLs export to match the SequenzD field +read orientation tag

### 0.2.3

Confindr "Inter" and "Inter+Intra" contamination now only result in a warning. This warning will transform in a fail only in combination with SNV or kraken fails.
THis behavior is consistent with the latest BfR thresholds.

### 0.2.2

Does not export fastas failing QC anymore

### 0.2.1

Relaxed regex for sample naming, subsample number now accepts tow to many digits (eg. 2023-0001254-01 or 2023-0001254-5236)

### 0.2.0

The workflow will expect fastq to be named after either the value of the `isolate_name_alt` if it is filled
or the value of the `isolate_id` field if `isolate_name_alt` is empty. If a fastq has a name that is not in these fields,
the workflow will stop with an error.

### 0.1.1

Add missing metadata field "assembly_method" in metadata JSON for geuebt export

### 0.1.0

First release

