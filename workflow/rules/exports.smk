# Read in data form aquamis summary and create data nescessary
# for export to geuebt, NRLs, etc...


rule geuebt_assemblies:
    input:
        ass="aquamis/Assembly/assembly/{sample_id}.fasta",
    output:
        ass="geuebt_export/{sample_id}.fasta",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/geuebt_assemblies_{sample_id}.log",
    shell:
        """
        exec 2> {log}
        cp {input.ass} {output.ass}
        """


rule geuebt_metadata:
    input:
        assembly=aggregate_assemblies,
        summary="aquamis/reports/summary_report.tsv",
    output:
        metatbl="geuebt_export/metadata.tsv",
    params:
        metadata=config["metadata"],
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/geuebt_metadata.log",
    script:
        "../scripts/geuebt_metadata.py"


rule nrls_export:
    input:
        # using geuebt table here becaus it's already qc checked
        metadata="geuebt_export/metadata.tsv",
        ssheet="sample_sheet/samples.tsv",
    output:
        outdir=directory("nrls_export"),
        flag=touch("nrls_export/sucess.flag"),
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/nrls_export.log",
    script:
        "../scripts/nrls_export.py"
