import os
import sys

from tempfile import TemporaryDirectory
import shutil
import filecmp
from pathlib import Path, PurePosixPath


sys.path.insert(0, os.path.dirname(__file__))


def test_geuebt_export():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/geuebt_export/data")
        expected_path = PurePosixPath(".tests/unit/geuebt_export/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/geuebt_export.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)
        os.mkdir(os.path.join(Path(tmpdir), "fasta_dump"))

        # run function
        sys.path.insert(0, workdir)
        from geuebt_export import main  # import main from your script
        main(
            summary=os.path.join(workdir, "summary_report.tsv"),
            metadata_in=os.path.join(workdir, "metadata.tsv"),
            assembly_path=workdir,
            fasta_dest=os.path.join(Path(tmpdir), "fasta_dump"),
            metadata_out=os.path.join(workdir, 'result.tsv'),
        )

        # check that tables are same
        assert filecmp.cmp(os.path.join(workdir, 'result.tsv'), os.path.join(expected_path, 'metadata.tsv'))
