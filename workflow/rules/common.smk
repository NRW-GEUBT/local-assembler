import os
import time
import pandas as pd
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


def validate_input_param(path, schema):
    df = pd.read_csv(path, index_col=False, sep="\t", engine="python")
    validate(df, schema=schema)


# Input functions ------------------------------------
def aggregate_assemblies(wildcards):
    checkpoint_output = checkpoints.aquamis.get(**wildcards).output["assdir"]
    ids_map = glob_wildcards(
        os.path.join(checkpoint_output, "{sample_id}.fasta")
    ).sample_id
    return expand("geuebt_export/{sample_id}.fasta", sample_id=ids_map)


# Validation
validate(config, schema="../schema/config.schema.yaml")

