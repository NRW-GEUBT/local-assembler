import os
import sys

from tempfile import TemporaryDirectory
import shutil
import filecmp
from pathlib import Path, PurePosixPath


sys.path.insert(0, os.path.dirname(__file__))


def test_geuebt_metadata():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/geuebt_metadata/data")
        expected_path = PurePosixPath(".tests/unit/geuebt_metadata/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/geuebt_metadata.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from geuebt_metadata import main  # import main from your script
        main(
            assemblies=[os.path.join(workdir, "2016-0000962-01.fasta")],
            summary=os.path.join(workdir, "summary_report.tsv"),
            metadata_in=os.path.join(workdir, "metadata.tsv"),
            metadata_out=os.path.join(workdir, 'result.tsv'),
        )

        # check that tables are same
        assert filecmp.cmp(os.path.join(workdir, 'result.tsv'), os.path.join(expected_path, 'metadata.tsv'))
