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


rule consolidate_ids:
    input:
        sample_sheet="sample_sheet/samples.tsv",
    params:
        metadata=config["metadata"],
    output:
        sample_sheet="sample_sheet/samples_isolate_ids.tsv",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/consolidate_ids.log",
    script:
        "../scripts/consolidate_ids.py"


checkpoint aquamis:
    input:
        sample_sheet="sample_sheet/samples_isolate_ids.tsv",
    output:
        outdir=directory("aquamis"),
        summary="aquamis/reports/summary_report.tsv",
        assdir=directory("aquamis/Assembly/assembly"),
    params:
        max_threads_sample=config["max_threads_sample"],
        qc_schema=f"{workflow.basedir}/schema/AQUAMIS_thresholds.json",
        run_name=f"-r {config['run_name']}" if config["run_name"] else "",
    conda:
        "../envs/aquamis.yaml"
    threads: workflow.cores
    log:
        "logs/run_aquamis.log",
    shell:
        """
        exec 2> {log}
        aquamis --sample_list {input.sample_sheet} \
            --working_directory {output.outdir} \
            --threads {threads} \
            --threads_sample {params.max_threads_sample} \
            --remove_temp --ephemeral \
            --qc_thresholds {params.qc_schema} \
            --fix_fails {params.run_name}
        """
