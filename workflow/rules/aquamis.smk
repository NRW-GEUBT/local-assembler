# Create sample sheet, run aquamis and validate qc values


rule create_sample_sheet:
    output:
        outdir=directory("sample_sheet"),
        sample_sheet="sample_sheet/samples.tsv",
    params:
        fastq_folder=config["fastq_folder"],
        fastq_naming=config["fastq_naming"],
    conda:
        "../envs/aquamis.yaml"
    log:
        "logs/create_sample_sheet.log",
    shell:
        """
        exec 2> {log}
        create_sampleSheet.sh \
            --mode {params.fastq_naming} \
            --fastxDir $(realpath {params.fastq_folder}) \
            --outDir {output.outdir} \
            --force
        """


checkpoint aquamis:
    input:
        sample_sheet="sample_sheet/samples.tsv",
    output:
        outdir=directory("aquamis"),
        summary="aquamis/reports/summary_report.tsv",
        assdir=directory("aquamis/Assembly/assembly"),
    params:
        max_threads_sample=config["max_threads_sample"],
        qc_schema=f"{workflow.basedir}/schema/AQUAMIS_thresholds.json",
    conda:
        "../envs/aquamis.yaml"
    threads: workflow.cores
    log:
        "logs/run_aquamis.log",
    shell:
        """
        aquamis --sample_list {input.sample_sheet} \
            --working_directory {output.outdir} \
            --threads {threads} \
            --threads_sample {params.max_threads_sample} \
            --remove_temp \
            --qc_thresholds {params.qc_schema}
        """
