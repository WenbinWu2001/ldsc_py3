"""Import-direction rules for the ``ldsc`` package layout.

The package is split into private domain packages arranged in tiers. A module
may import from its own package and from lower tiers only, and the peer domains
in one tier may not import each other. The check reads every import inside the
package, including imports inside functions and under ``TYPE_CHECKING``, so a
function-level import cannot hide an upward dependency. Every module outside
the root shell (``__init__``, ``__main__``, ``cli``, ``config`` and ``errors``)
belongs to a tiered package.
"""

import ast
import unittest
from pathlib import Path

import ldsc

# Root modules are matched exactly; packages also own every module below them.
SHELL_TIERS = {
    "ldsc.errors": 0,
    "ldsc.config": 3,
    "ldsc.cli": 9,
    "ldsc.__main__": 9,
    "ldsc": 9,
}
# The entry points in the top tier may import any module.
TOP_TIER = 9
PACKAGE_TIERS = {
    "ldsc._core": 1,
    "ldsc._runtime": 2,
    "ldsc._genome": 4,
    "ldsc._annotation": 5,
    "ldsc._refpanel": 5,
    "ldsc._sumstats": 5,
    "ldsc._ldscore": 6,
    "ldsc._regression": 7,
    "ldsc._plotting": 8,
}


def _module_name(package_root, path):
    parts = list(path.relative_to(package_root.parent).with_suffix("").parts)
    return ".".join(parts[:-1] if parts[-1] == "__init__" else parts)


def _unit(module):
    """Return the tiered unit that owns ``module``, or ``None`` when it has no tier yet."""
    if module in SHELL_TIERS:
        return module
    for package in PACKAGE_TIERS:
        if module == package or module.startswith(package + "."):
            return package
    return None


def _tier(unit):
    return SHELL_TIERS.get(unit, PACKAGE_TIERS.get(unit))


def _imported_modules(tree, module, is_package, known):
    package = module if is_package else module.rpartition(".")[0]
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "ldsc" or alias.name.startswith("ldsc."):
                    yield node.lineno, alias.name
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.split(".")
                base = base[: len(base) - (node.level - 1)]
                target = ".".join(base + ([node.module] if node.module else []))
            elif node.module and (node.module == "ldsc" or node.module.startswith("ldsc.")):
                target = node.module
            else:
                continue
            for alias in node.names:
                submodule = f"{target}.{alias.name}"
                yield node.lineno, submodule if submodule in known else target


def find_violations(package_root):
    """Return every import in ``package_root`` that breaks the tier rules."""
    package_root = Path(package_root)
    paths = sorted(package_root.rglob("*.py"))
    known = {_module_name(package_root, path) for path in paths}
    violations = []
    for path in paths:
        module = _module_name(package_root, path)
        unit = _unit(module)
        if unit is None or _tier(unit) == TOP_TIER:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports = set(_imported_modules(tree, module, path.name == "__init__.py", known))
        for line, target in sorted(imports):
            target_unit = _unit(target)
            where = f"{path.relative_to(package_root.parent)}:{line} imports {target}"
            if target_unit is None:
                if unit in PACKAGE_TIERS:
                    violations.append(f"{where}, which is not yet in a tiered package")
            elif target_unit != unit:
                if _tier(target_unit) > _tier(unit):
                    violations.append(f"{where} (tier {_tier(target_unit)}) from tier {_tier(unit)}")
                elif _tier(target_unit) == _tier(unit):
                    violations.append(f"{where}, a peer domain in tier {_tier(unit)}")
    return violations


def find_misplaced_modules(package_root):
    """Return modules and packages at the package root outside the shell and the tiered packages."""
    package_root = Path(package_root)
    shell = {name.rpartition(".")[2] for name in SHELL_TIERS if name != "ldsc"} | {"__init__"}
    packages = {name.rpartition(".")[2] for name in PACKAGE_TIERS}
    misplaced = []
    for path in sorted(package_root.iterdir()):
        if path.suffix == ".py" and path.stem not in shell:
            misplaced.append(path.name)
        elif path.is_dir() and path.name not in packages and any(path.rglob("*.py")):
            misplaced.append(path.name + "/")
    return misplaced


class LayeringTest(unittest.TestCase):
    def test_imports_point_to_lower_tiers_only(self):
        violations = find_violations(Path(ldsc.__file__).parent)
        self.assertEqual(violations, [], "\n" + "\n".join(violations))

    def test_package_root_holds_only_the_shell(self):
        misplaced = find_misplaced_modules(Path(ldsc.__file__).parent)
        self.assertEqual(misplaced, [], "Move these into a tiered package: " + ", ".join(misplaced))


if __name__ == "__main__":
    unittest.main()
