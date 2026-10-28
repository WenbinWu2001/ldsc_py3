from __future__ import annotations

import subprocess
import sys

def test_root_import_and_lazy_public_exports_do_not_import_matplotlib():
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys; import ldsc; "
                "assert 'matplotlib' not in sys.modules; "
                "from ldsc import plot_result, convert_h2_scale; "
                "assert 'matplotlib' not in sys.modules"
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
