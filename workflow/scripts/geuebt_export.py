#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import hashlib
import pandas as pd


def md5(fname):
    """
    Return the hex string representation of the md5 difest for a file
    From stackoverflow user @quantumSoup (2010-08-7)
    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def main(summary, metadata_in, assembly_path, fasta_dest, metadata_out):
    qc = []
    # load metadata and summary
    metatbl = pd.read_csv(metadata_in, sep="\t", index_col="isolate_id")
    sumtbl = pd.read_csv(summary, sep="\t", index_col="Sample_Name")
    # Insert assembly method
    metatbl['assembly_method'] = "AQUAMIS"
    # for each assembly
    for name in sumtbl.index:
        fastaname = os.path.join(assembly_path, f"{name}.fasta")
        # Just crash the workflow with a message if metadata are missing
        if name not in metatbl.index:
            raise KeyError(
                f"There is not information on sample '{name}' in the metadata table, "
                f"althought valid FASTQs were provided. "
                f"Ensure the completness of the submitted metadata. "
                f"The workflow will expect fastq to be named after either the value "
                f"the `isolate_id` field if `isolate_name_alt` is empty."
            )
        # iof QC fail skip sample
        if sumtbl.at[name, "QC_Vote"] == "FAIL":
            continue
        # copy files
        os.symlink(fastaname, os.path.join(fasta_dest, f"{name}.fasta"))
        # QC values directly as a single row of a dataframe
        qc.append(pd.DataFrame.from_dict(
            {name: [
                os.path.basename(fastaname),
                md5(fastaname),
                sumtbl.at[name, "Coverage_Depth"],
                sumtbl.at[name, "Reference_Coverage"],
                sumtbl.at[name, "Q30_Base_Fraction"]
            ]},
            orient="index",
            columns=[
                "fasta_name",
                "fasta_md5",
                "seq_depth",
                "ref_coverage",
                "q30"
            ]
        ))
    # dict to df
    metaext = pd.concat(qc)
    # merge tables (inner join)
    merged = pd.merge(metatbl, metaext, how="inner", left_index=True, right_index=True)
    merged.index.name = "isolate_id"
    # output
    merged.to_csv(metadata_out, sep="\t", header=True, index=True)


if __name__ == "__main__":
    main(
        snakemake.input["summary"],
        snakemake.params["metadata"],
        snakemake.params["assembly_path"],
        snakemake.params["fasta_destination"],
        snakemake.output["metatbl"],
    )
