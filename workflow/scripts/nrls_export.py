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
import shutil
import datetime
import pandas as pd


EXPORT_CONDITIONS = ["lebensmittel"]


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


def main(metadata, ssheet, outdir):
    # make outdir
    os.makedirs(outdir, exist_ok=True)
    # load metadata and ssheet as df
    metatbl = pd.read_csv(metadata, sep="\t", index_col="isolate_id")
    # Only for subset of sample types !!!
    selected_meta = metatbl
    # selected_meta = metatbl[metatbl["sample_type"].isin(EXPORT_CONDITIONS)]
    ssheet = pd.read_csv(ssheet, sep="\t", index_col="sample")
    # left join on metadata
    tbl = pd.merge(selected_meta, ssheet, how="left", left_index=True, right_index=True)
    # for each species:
    for species in tbl["organism"].unique():
        checksums, metanrl = [], []
        species = species.replace(" ", "_").replace(".", "")
        # create a subfolder
        os.makedirs(os.path.join(outdir, species), exist_ok=True)
        # for each sample:
        for row in tbl.iterrows():  # yields (index, Series)
            for fastqpath in (row[1]["fq1"], row[1]["fq2"]):
                filename = os.path.basename(fastqpath)
                # rename files with isolate_id
                filename_id = filename.replace(row[1]["fastq_name"], row[0])
                # copy fastq
                shutil.copy(fastqpath, os.path.join(outdir, species, filename_id))
                # get checksum
                checksums.append(f"{md5(fastqpath)}  {filename_id}")
            # create a one row df for metadata
            metanrl.append(pd.DataFrame.from_dict(
                {row[0]: [
                    row[0],  # "SequenzID_Einsender",
                    "fastq",  # "Art_Sequenzdaten",
                    row[1]["sequencing_instrument"],  # "Sequenziergerät",
                    "<FILL IN CREDENTIALS>",  # "EinsenderID",
                    "<FILL IN CREDENTIALS>",  # "Einsender_Institution",
                    "<FILL IN CREDENTIALS>",  # "Einsender_Email_Ansprechperson",
                    datetime.date.today(),  # "Datum_SequenzUpload",
                    "",  # "Zweck_Datenaustausch",
                    "",  # "Isolatversand_an_BfR_(ja/nein)",
                    "",  # "Datum_Isolatversand",
                    "",  # "BfR_Probennummer(falls bekannt)",
                    row[1]["sample_id"],  # "Original-Nr._Probe",
                    row[1]["sample_id"],  # "AVVData-Nr._Probe",
                    species,  # "Erreger",
                    "",  # "Serovar",
                    row[1]["collection_date"],  # "Probenahme_Datum",
                    "",  # "Isolierung_Datum",
                    "",  # "Probenahme_Ort_Code(ADV-Kat-Nr.9)",
                    "",  # "Probenahme_Ort_Text(ADV-Kat-Nr.9)",
                    "",  # "Probenahme_Ort(PLZ)",
                    row[1]["collection_municipality"],  # "Probenahme_Ort(STADT)",
                    "",  # "Matrix_Oberbegriff_Code(ADV-Kat-Nr.2)",
                    "",  # "Matrix_Code(ADV-Kat-Nr.3)",
                    row[1]["designation"],  # "Matrix_Text(ADV-Kat-Nr.3)",
                    "",  # "Matrix_Textfeld_Ergänzung",
                    "",  # "Probenahme_Grund_Code(ADV-Kat-Nr.4)",
                    row[1]["collection_cause"],  # "Probenahme_Grund_Text(ADV-Kat-Nr.4)",
                    "",  # "Probenahme_Grund_Textfeld_Ergänzung",
                    "",  # "Betriebsart_Code(ADV-Kat-Nr.8)",
                    "",  # "Betriebsart_Text(ADV-Kat-Nr.8)",
                    "",  # "Betriebsart_Textfeld_Ergänzung",
                    "",  # "Probe_Bemerkung",
                ]},
                orient="index",
                columns=[
                    "SequenzID_Einsender",
                    "Art_Sequenzdaten",
                    "Sequenziergerät",
                    "EinsenderID",
                    "Einsender_Institution",
                    "Einsender_Email_Ansprechperson",
                    "Datum_SequenzUpload",
                    "Zweck_Datenaustausch",
                    "Isolatversand_an_BfR_(ja/nein)",
                    "Datum_Isolatversand",
                    "BfR_Probennummer(falls bekannt)",
                    "Original-Nr._Probe",
                    "AVVData-Nr._Probe",
                    "Erreger",
                    "Serovar",
                    "Probenahme_Datum",
                    "Isolierung_Datum",
                    "Probenahme_Ort_Code(ADV-Kat-Nr.9)",
                    "Probenahme_Ort_Text(ADV-Kat-Nr.9)",
                    "Probenahme_Ort(PLZ)",
                    "Probenahme_Ort(STADT)",
                    "Matrix_Oberbegriff_Code(ADV-Kat-Nr.2)",
                    "Matrix_Code(ADV-Kat-Nr.3)",
                    "Matrix_Text(ADV-Kat-Nr.3)",
                    "Matrix_Textfeld_Ergänzung",
                    "Probenahme_Grund_Code(ADV-Kat-Nr.4)",
                    "Probenahme_Grund_Text(ADV-Kat-Nr.4)",
                    "Probenahme_Grund_Textfeld_Ergänzung",
                    "Betriebsart_Code(ADV-Kat-Nr.8)",
                    "Betriebsart_Text(ADV-Kat-Nr.8)",
                    "Betriebsart_Textfeld_Ergänzung",
                    "Probe_Bemerkung",
                ]
            ))
        # concatenate df
        metaext = pd.concat(metanrl)
        # output both tables
        with open(os.path.join(outdir, species, f"md5_{species}.txt"), "w") as fo:
            fo.write("\n".join(checksums))
        metaext.to_csv(
            os.path.join(outdir, species, f"metadata_{species}.tsv"),
            sep="\t",
            header=True,
            index=False
        )


if __name__ == "__main__":
    main(
        snakemake.input["metadata"],
        snakemake.input["ssheet"],
        snakemake.output["outdir"],
    )
