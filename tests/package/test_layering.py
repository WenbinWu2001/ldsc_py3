"""Import-direction rules for the ``ldsc`` package layout.

The package is split into private domain packages arranged in tiers. A module
may import from its own package and from lower tiers only, and the peer domains
in one tier may not import each other. The check reads every import inside the
package, including imports inside functions and under ``TYPE_CHECKING``, so a
function-level import cannot hide an upward dependency. Every module outside
the root shell (``__init__``, ``__main__``, ``cli``, ``config`` and ``errors``)
belongs to a tiered package.

Two further rules govern the import graph Python actually runs, which counts
function-level imports and leaves out ``TYPE_CHECKING`` imports: no modules may
import each other, directly or through others, and an in-package import may sit
inside a function only when ``DEFERRED_IMPORTS`` lists it with its reason.
"""

import ast
import json
import subprocess
import sys
import tempfile
import textwrap
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

# The only in-package imports that stay inside functions, as
# (importing module, imported module): reason. `import ldsc` or `ldsc --help`
# loads each importer but not its target (checked below), so moving the
# import to module level would load the target at every start. The Matplotlib
# plot builders follow their own rule instead: Matplotlib loads only on paths
# that draw a figure, so `ldsc plot --help`, which loads their importer, must
# not load them either (FIGURE_IMPORTS, also checked below).
FIGURE_IMPORTS = {("ldsc._plotting.dispatch", "ldsc._plotting.builders")}
DEFERRED_IMPORTS = {
    ("ldsc._plotting.dispatch", "ldsc._plotting.builders"):
        "the plot builders import Matplotlib, which loads only for a figure",
    ("ldsc._annotation.bundle", "ldsc._annotation.in_memory"):
        "only AnnotationBundle.from_frames builds in-memory annotations",
    ("ldsc._ldscore.calculator", "ldsc._annotation.in_memory"):
        "only in-memory Python runs keep dropped-SNP diagnostics in memory",
    ("ldsc._ldscore.calculator", "ldsc._ldscore.direct_annotation"):
        "only the ldscore command prepares direct BED and gene-list annotations",
    ("ldsc._ldscore.calculator", "ldsc._ldscore.overlap_table"):
        "the overlap matrix is assembled only after chromosome results exist",
    ("ldsc._ldscore.calculator", "ldsc._ldscore.preflight"):
        "only the ldscore command preflights its reference inputs",
    ("ldsc._ldscore.calculator", "ldsc._refpanel.preflight"):
        "only R2 matrix runs check their sidecar headers",
    ("ldsc._ldscore.calculator", "ldsc._runtime.workflow_kwargs"):
        "only the ldscore entry points check keyword arguments and record explicit options",
    ("ldsc._ldscore.gene_index.build", "ldsc._refpanel.preflight"):
        "only an index build checks genetic-map headers",
    ("ldsc._ldscore.gene_index.scoring", "ldsc._ldscore.gene_index.batches"):
        "only indexed scoring runs the indexed LD-score producer",
    ("ldsc._ldscore.gene_index.scoring", "ldsc._ldscore.preflight"):
        "only a failed index load writes an input-issue table",
    ("ldsc._ldscore.outputs", "ldsc._ldscore.overlap_table"):
        "only directories with an overlap matrix write one",
    ("ldsc._ldscore.source", "ldsc._ldscore.overlap_table"):
        "only directories with an overlap matrix read one",
    ("ldsc._refpanel.builder", "ldsc._refpanel.preflight"):
        "only build-r2-matrix checks genetic-map headers",
    ("ldsc._refpanel.builder", "ldsc._runtime.workflow_kwargs"):
        "only run_build_r2_matrix checks Python keyword arguments",
    ("ldsc._refpanel.query", "ldsc._refpanel.preflight"):
        "only query-r2-matrix checks R2 matrix sidecar headers",
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


def _is_type_checking(test):
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING")


def _runtime_imports(node, in_function=False):
    """Yield ``(import node, inside a function)`` for imports outside ``TYPE_CHECKING`` blocks."""
    for child in ast.iter_child_nodes(node):
        if isinstance(child, ast.If) and _is_type_checking(child.test):
            for item in child.orelse:
                yield from _runtime_imports_of(item, in_function)
            continue
        yield from _runtime_imports_of(child, in_function)


def _runtime_imports_of(node, in_function):
    if isinstance(node, (ast.Import, ast.ImportFrom)):
        yield node, in_function
    yield from _runtime_imports(node, in_function or isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)))


def _runtime_edges(package_root):
    """Return ``(module, line, imported module, inside a function)`` for every runtime in-package import."""
    package_root = Path(package_root)
    paths = sorted(package_root.rglob("*.py"))
    known = {_module_name(package_root, path) for path in paths}
    edges = []
    for path in paths:
        module = _module_name(package_root, path)
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node, in_function in _runtime_imports(tree):
            for line, target in _imported_modules(node, module, path.name == "__init__.py", known):
                if target in known:
                    edges.append((module, line, target, in_function))
    return edges


def find_import_cycles(package_root):
    """Return each group of modules that import one another at run time, as sorted lists.

    Function-level imports count; ``TYPE_CHECKING`` imports do not.
    """
    graph = {}
    for module, _, target, _ in _runtime_edges(package_root):
        graph.setdefault(module, set()).add(target)
        graph.setdefault(target, set())
    # Tarjan's strongly connected components, iteratively.
    index, low, on_stack, stack, cycles = {}, {}, set(), [], []
    for start in sorted(graph):
        if start in index:
            continue
        work = [(start, iter(sorted(graph[start])))]
        index[start] = low[start] = len(index)
        stack.append(start)
        on_stack.add(start)
        while work:
            node, children = work[-1]
            child = next(children, None)
            if child is None:
                work.pop()
                if work:
                    low[work[-1][0]] = min(low[work[-1][0]], low[node])
                if low[node] == index[node]:
                    component = []
                    while True:
                        member = stack.pop()
                        on_stack.discard(member)
                        component.append(member)
                        if member == node:
                            break
                    if len(component) > 1 or node in graph[node]:
                        cycles.append(sorted(component))
            elif child not in index:
                index[child] = low[child] = len(index)
                stack.append(child)
                on_stack.add(child)
                work.append((child, iter(sorted(graph[child]))))
            elif child in on_stack:
                low[node] = min(low[node], index[child])
    return sorted(cycles)


def find_function_level_imports(package_root):
    """Map each runtime in-package import inside a function to the places it occurs."""
    package_root = Path(package_root)
    found = {}
    for module, line, target, in_function in _runtime_edges(package_root):
        if in_function:
            found.setdefault((module, target), []).append(f"{module}:{line}")
    return found


def find_misplaced_modules(package_root):
    """Return modules and packages at the package root outside the root shell and the tiered packages."""
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

    def test_no_runtime_import_cycles(self):
        cycles = find_import_cycles(Path(ldsc.__file__).parent)
        self.assertEqual(cycles, [], "These modules import one another; move the shared code so imports point "
                                     "one way:\n" + "\n".join(", ".join(cycle) for cycle in cycles))

    def test_function_level_imports_are_listed_with_a_reason(self):
        found = find_function_level_imports(Path(ldsc.__file__).parent)
        unlisted = sorted(place for pair, places in found.items() if pair not in DEFERRED_IMPORTS for place in places)
        self.assertEqual(unlisted, [], "Move these in-package imports to module level, or add the pair to "
                                       "DEFERRED_IMPORTS with its reason:\n" + "\n".join(unlisted))
        stale = sorted(f"{importer} -> {target}" for importer, target in DEFERRED_IMPORTS
                       if (importer, target) not in found)
        self.assertEqual(stale, [], "DEFERRED_IMPORTS lists imports that no longer exist:\n" + "\n".join(stale))
        self.assertTrue(all(reason.strip() for reason in DEFERRED_IMPORTS.values()))

    def test_deferred_imports_stay_unloaded_at_startup(self):
        # A clean interpreter shows whether the premise of each deferred import
        # still holds: at startup its importer is loaded and its target is not,
        # and the Matplotlib builders stay unloaded by `ldsc plot --help` too.
        probe = textwrap.dedent(
            """
            import contextlib, io, json, sys
            import ldsc
            from ldsc import cli

            def loaded_after(argv):
                with contextlib.redirect_stdout(io.StringIO()):
                    try:
                        cli.main(argv)
                    except SystemExit:
                        pass
                return sorted(name for name in sys.modules if name.startswith("ldsc"))

            print(json.dumps({"startup": loaded_after(["--help"]), "plot_help": loaded_after(["plot", "--help"])}))
            """
        )
        completed = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True,
                                   cwd=tempfile.gettempdir())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        loaded = {name: set(modules) for name, modules in json.loads(completed.stdout.splitlines()[-1]).items()}
        startup_targets = sorted(target for _, target in DEFERRED_IMPORTS if target in loaded["startup"])
        self.assertEqual(startup_targets, [],
                         "`import ldsc` or `ldsc --help` loads a deferred module, so its imports no longer "
                         "need to be deferred")
        idle_importers = sorted(f"{importer} -> {target}" for importer, target in DEFERRED_IMPORTS
                                if (importer, target) not in FIGURE_IMPORTS and importer not in loaded["startup"])
        self.assertEqual(idle_importers, [],
                         "`import ldsc` and `ldsc --help` do not load these importers, so a module-level "
                         "import would not slow startup; move it to module level:\n" + "\n".join(idle_importers))
        for importer, target in FIGURE_IMPORTS:
            self.assertIn(importer, loaded["plot_help"])
            self.assertNotIn(target, loaded["plot_help"], "`ldsc plot --help` loads the Matplotlib plot builders")

if __name__ == "__main__":
    unittest.main()
