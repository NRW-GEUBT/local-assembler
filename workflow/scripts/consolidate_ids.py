#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import pandas as pd


def main(ssheet, metadata, sheetout):
    # load metadata and samplesheet
    metatbl = pd.read_csv(metadata, sep="\t", index_col=False)
    sample_sheet = pd.read_csv(ssheet, sep="\t", index_col="sample")
    new_index = []
    # for each fastq pair, check that there is a corresponding entry in metadata
    for sname in sample_sheet.index.to_list():
        selection = metatbl.loc[metatbl["alt_isolate_id"] == sname]
        if len(selection) == 0:
            selection = metatbl.loc[metatbl["isolate_id"] == sname]
        if len(selection) == 0:
            # Crash if name not found
            raise KeyError(
                f"There is not information on sample '{sname}' in the metadata table, "
                f"althought valid FASTQs were provided. "
                f"Ensure the completness of the submitted metadata. "
                f"The workflow will expect fastq to be named after either the value "
                f"the `isolate_id` field if `alt_isolate_id` is empty."
            )
        elif len(selection) > 1:
            # crash if name not unique
            raise ValueError(
                f"Several entries for sample name '{sname}' were found in the metadata table. "
                f" Make sure that both the fields `isolate_id` and `alt_isolate_id` "
                f"contain unique values"
            )
        else:
            # add isolate id to reindexing
            new_index.append(selection.iloc[0]["isolate_id"])

    # resindex
    sample_sheet.reset_index(inplace=True, names="fastq_name")
    reindexed_sheet = sample_sheet.set_index(pd.Index(new_index))
    reindexed_sheet.index.name = "sample"
    # output
    reindexed_sheet.to_csv(sheetout, sep="\t", header=True, index=True)


if __name__ == "__main__":
    main(
        snakemake.input["sample_sheet"],
        snakemake.params["metadata"],
        snakemake.output["sample_sheet"],
    )
