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

