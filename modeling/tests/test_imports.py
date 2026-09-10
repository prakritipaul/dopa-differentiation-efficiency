"""Every module must import. Cheap guard for entry points no other test touches.

modeling/run_variant.py once shipped with an IndentationError: it is an
entry point, no test imported it, so nothing noticed until it was run by
hand. Importing catches syntax errors, bad imports, and anything that
breaks at module scope -- which is where path constants are defined.
"""

import importlib
import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
MODULES = sorted(p.stem for p in (REPO_ROOT / "modeling").glob("*.py") if p.stem != "__init__")
SCRIPTS = sorted(p.name for p in REPO_ROOT.glob("0*.py"))


@pytest.mark.parametrize("name", MODULES)
def test_modeling_module_imports(name):
    importlib.import_module(f"modeling.{name}")


@pytest.mark.parametrize("name", SCRIPTS)
def test_numbered_script_imports(name):
    # Numbered EDA scripts aren't a package; load them by path the same way
    # modeling/features.py does. main() is __main__-guarded, so importing
    # only executes module-level definitions.
    spec = importlib.util.spec_from_file_location(name[:-3], REPO_ROOT / name)
    spec.loader.exec_module(importlib.util.module_from_spec(spec))


def test_numbered_script_output_paths_resolve():
    """Every OUT_DIR / "..." literal in the EDA scripts must land in a real
    directory.

    These scripts are never run by the test suite -- they stream multi-GB h5
    files -- so a reorganisation can silently repoint their writes. Checking
    the literals statically is cheap and catches the whole failure mode:
    a wrong prefix, or a missing one that would dump output at the
    metadata_eda/ root instead of its subdirectory.
    """
    import re

    out_dir = REPO_ROOT / "metadata_eda"
    broken = []
    for script in REPO_ROOT.glob("0*.py"):
        for m in re.finditer(r'OUT_DIR\s*/\s*f?"([^"]+\.(?:csv|png))"', script.read_text()):
            literal = m.group(1)
            if "/" not in literal:
                broken.append(f"{script.name}: {literal!r} has no subdirectory prefix")
            elif not (out_dir / literal).parent.is_dir():
                broken.append(f"{script.name}: {literal!r} -> parent directory does not exist")
    assert not broken, "output paths would not resolve:\n  " + "\n  ".join(broken)
