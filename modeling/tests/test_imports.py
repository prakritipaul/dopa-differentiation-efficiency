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
