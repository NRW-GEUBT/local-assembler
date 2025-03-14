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
import datetime
import pandas as pd


# Only these will be exported, third party always excluded
EXPORT_CONDITIONS = ["Lebensmittel", "Umfeld", "Futtermittel", "Tiergesundheit"]


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
    selected_meta = metatbl.loc[(metatbl["sample_type"].isin(EXPORT_CONDITIONS)) & (metatbl["third_party_owner"]=="")]
    ssheet = pd.read_csv(ssheet, sep="\t", index_col="sample")
    # left join on metadata
    tbl = pd.merge(selected_meta, ssheet, how="left", left_index=True, right_index=True)
    # for each species:
    for species in tbl["organism"].unique():
        checksums, metanrl = [], []
        species_fmt = species.replace(" ", "_").replace(".", "")
        # create a subfolder
        os.makedirs(os.path.join(outdir, species_fmt), exist_ok=True)
        # select by org
        subtbl = tbl.loc[tbl["organism"] == species]
        # for each sample:
        for row in subtbl.iterrows():  # yields (index, Series)
            for fastqpath in zip([row[1]["fq1"], row[1]["fq2"]], ["R1", "R2"]):
                # rename files with isolate_id
                renamed = os.path.join(outdir, species_fmt, f"{row[0]}_{fastqpath[1]}.fastq.gz")

                # Symlink fastq
                try:
                    os.symlink(os.path.abspath(fastqpath[0]), renamed)
                except FileExistsError:
                    os.remove(renamed)
                    os.symlink(os.path.abspath(fastqpath[0]), renamed)

                # get checksum
                checksums.append(f"{md5(renamed)}  {row[0]}_{fastqpath[1]}.fastq.gz")

            # create a one row df for metadata
            metanrl.append(pd.DataFrame.from_dict(
                {row[0]: [
                    row[0], # SequenzID_Einsender
                    "fastq", # Art_Sequenzdaten
                    row[1]["sequencing_instrument"], # Sequenziergerät
                    "<FILL IN ID>", # EinsenderID
                    "<FILL IN INSTITUTION>", # Einsender_Institution
                    "<FILL IN EMAIL>", # Einsender_Email_Ansprechperson
                    datetime.date.today(), # Datum_SequenzUpload
                    "Abgleich", # Zweck_Datenaustausch
                    "ja", # Isolatversand_an_BfR_(ja/nein)
                    "", # Datum_Isolatversand
                    "", # BfR-ProbenNr/falls_bekannt)
                    row[1]["sample_id"], # Orig.-Nr.
                    row[1]["sample_id"], # AVVDatA-ProbenNr
                    row[1]["isolate_id"], # AVVDatA-TeilprobenNr
                    row[1]["organism"], # Erreger
                    "", # Serovar
                    row[1]["collection_date"], # Probenahmedatum
                    row[1]["isolation_date"], # Isolierdatum
                    row[1]["collection_cause_code"], # Probenahmegrund(Code_aus_ADV-Kat-Nr.4)
                    row[1]["collection_cause"], # Probenahmegrund(Text_aus_ADV-Kat-Nr.4)
                    row[1]["control_program_code"], # Kontrollprogramm(Code_aus_AVVDatA-Kat-Nr.322)
                    row[1]["control_program"], # Kontrollprogramm(Text_aus_AVVDatA-Kat-Nr.322)
                    row[1]["analysis_cause_code"], # Untersuchungsgrund(Code_aus_AVVDatA-Kat-Nr.326)
                    row[1]["analysis_cause"], # Untersuchungsgrund(Text_aus_AVVDatA-Kat-Nr.326)
                    row[1]["control_program_description"], # Kontrollprogramm/Untersuchungsgrund(Freitext)
                    row[1]["matrix_group_code"], # Oberbegriff_Kodiersystem_der_Matrizes(Code_aus_ADV-Kat-Nr.2)
                    row[1]["matrix_code"], # Matrix(Code_aus_ADV-Kat-Nr.3)
                    row[1]["matrix"], # Matrix(Text_aus_ADV-Kat-Nr.3)
                    "", # Tiere(Code_aus_AVVDatA-Kat-Nr.339)
                    "", # Tiere(Text_aus_AVVDatA-Kat-Nr.339)
                    row[1]["matrix_code"], # Matrix(Code_aus_AVVDatA-Kat-Nr.319)
                    row[1]["matrix"], # Matrix(Text_aus_AVVDatA-Kat-Nr.319)
                    row[1]["description"], # Tiere/Matrix(Freitext)
                    "", # Angaben_zur_Primärproduktion(Code_aus_AVVDatA-Kat-Nr.316)
                    "", # Angaben_zur_Primärproduktion(Text_aus_AVVDatA-Kat-Nr.316)
                    row[1]["collection_place_code"], # Ort_der_Probenahme(Code_aus_ADV-Kat-Nr.9)
                    row[1]["collection_place"], # Ort_der_Probenahme(Text_aus_ADV-Kat-Nr.9)
                    row[1]["collection_place_code"], # Ort_der_Probenahme(Code_aus_AVVDatA-Kat-Nr.313)
                    row[1]["collection_place"], # Ort_der_Probenahme(Text_aus_AVVDatA-Kat-Nr.313)
                    row[1]["collection_place_city"], # Probenahmeort_STADT(Freitext)
                    row[1]["collection_place_zip"], # Probenahmeort_PLZ(Freitext)
                    row[1]["manufacturer_type_code"], # Betriebsart(Code_aus_ADV-Kat-Nr.8)
                    row[1]["manufacturer_type"], # Betriebsart(Text_aus_ADV-Kat-Nr.8)
                    row[1]["manufacturer_type_code"], # Betriebsart(Code_aus_AVVDatA-Kat-Nr.303)
                    row[1]["manufacturer_type"], # Betriebsart(Text_aus_AVVDatA-Kat-Nr.303)
                    row[1]["manufacturer_type_description"], # Betriebsart(Freitext)
                    "", # Programm(Code_aus_AVVDatA-Kat-Nr.328)
                    "", # Programm(Text_aus_AVVDatA-Kat-Nr.328)
                    row[1]["notes"], # Bemerkung
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
                    "BfR-ProbenNr/falls_bekannt)",
                    "Orig.-Nr.",
                    "AVVDatA-ProbenNr",
                    "AVVDatA-TeilprobenNr",
                    "Erreger",
                    "Serovar",
                    "Probenahmedatum",
                    "Isolierdatum",
                    "Probenahmegrund(Code_aus_ADV-Kat-Nr.4)",
                    "Probenahmegrund(Text_aus_ADV-Kat-Nr.4)",
                    "Kontrollprogramm(Code_aus_AVVDatA-Kat-Nr.322)",
                    "Kontrollprogramm(Text_aus_AVVDatA-Kat-Nr.322)",
                    "Untersuchungsgrund(Code_aus_AVVDatA-Kat-Nr.326)",
                    "Untersuchungsgrund(Text_aus_AVVDatA-Kat-Nr.326)",
                    "Kontrollprogramm/Untersuchungsgrund(Freitext)",
                    "Oberbegriff_Kodiersystem_der_Matrizes(Code_aus_ADV-Kat-Nr.2)",
                    "Matrix(Code_aus_ADV-Kat-Nr.3)",
                    "Matrix(Text_aus_ADV-Kat-Nr.3)",
                    "Tiere(Code_aus_AVVDatA-Kat-Nr.339)",
                    "Tiere(Text_aus_AVVDatA-Kat-Nr.339)",
                    "Matrix(Code_aus_AVVDatA-Kat-Nr.319)",
                    "Matrix(Text_aus_AVVDatA-Kat-Nr.319)",
                    "Tiere/Matrix(Freitext)",
                    "Angaben_zur_Primärproduktion(Code_aus_AVVDatA-Kat-Nr.316)",
                    "Angaben_zur_Primärproduktion(Text_aus_AVVDatA-Kat-Nr.316)",
                    "Ort_der_Probenahme(Code_aus_ADV-Kat-Nr.9)",
                    "Ort_der_Probenahme(Text_aus_ADV-Kat-Nr.9)",
                    "Ort_der_Probenahme(Code_aus_AVVDatA-Kat-Nr.313)",
                    "Ort_der_Probenahme(Text_aus_AVVDatA-Kat-Nr.313)",
                    "Probenahmeort_STADT(Freitext)",
                    "Probenahmeort_PLZ(Freitext)",
                    "Betriebsart(Code_aus_ADV-Kat-Nr.8)",
                    "Betriebsart(Text_aus_ADV-Kat-Nr.8)",
                    "Betriebsart(Code_aus_AVVDatA-Kat-Nr.303)",
                    "Betriebsart(Text_aus_AVVDatA-Kat-Nr.303)",
                    "Betriebsart(Freitext)",
                    "Programm(Code_aus_AVVDatA-Kat-Nr.328)",
                    "Programm(Text_aus_AVVDatA-Kat-Nr.328)",
                    "Bemerkung",
                ]
            ))
        # concatenate df
        metaext = pd.concat(metanrl)
        # output both tables
        with open(os.path.join(outdir, species_fmt, f"md5_{species_fmt}.txt"), "w") as fo:
            fo.write("\n".join(checksums))
        metaext.to_csv(
            os.path.join(outdir, species_fmt, f"metadata_{species_fmt}.tsv"),
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
