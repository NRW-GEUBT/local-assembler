# Read in data form aquamis summary and create data nescessary
# for export to geuebt, NRLs, etc...


rule geuebt_export:
    input:
        summary="aquamis/reports/summary_report.tsv",
    output:
        metatbl="geuebt_export/metadata.tsv",
    params:
        metadata=config["metadata"],
        assembly_path="aquamis/Assembly/assembly/",
        fasta_destination=lambda w, output: os.path.split(output["metatbl"])[0],
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/geuebt_export.log",
    script:
        "../scripts/geuebt_export.py"


rule nrls_export:
    input:
        # using geuebt table here because it's already qc checked
        metadata="geuebt_export/metadata.tsv",
        ssheet="sample_sheet/samples_isolate_ids.tsv",
    output:
        outdir=directory("nrls_export"),
        flag=touch("nrls_export/sucess.flag"),
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/nrls_export.log",
    script:
        "../scripts/nrls_export.py"

