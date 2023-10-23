import os
import sys

from tempfile import TemporaryDirectory
import shutil
import filecmp
from pathlib import Path, PurePosixPath


sys.path.insert(0, os.path.dirname(__file__))


def test_consolidate_ids():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/consolidate_ids/data")
        expected_path = PurePosixPath(".tests/unit/consolidate_ids/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/consolidate_ids.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from consolidate_ids import main  # import main from your script
        main(
            ssheet=os.path.join(workdir, "sample_sheet.tsv"),
            metadata=os.path.join(workdir, "metadata.tsv"),
            sheetout=os.path.join(workdir, 'result.tsv'),
        )

        # check that tables are same
        assert filecmp.cmp(os.path.join(workdir, 'result.tsv'), os.path.join(expected_path, 'sample_sheet.tsv'))
