"""Check the installed packages against the conda environment file of the repository.

This runs before anything else when the mwf command is launched.
It is skipped in packaged installs (e.g. bioconda) since they do not include the environment file.
It only uses the standard library and yaml so it does not depend on the packages it checks.
"""
import re
import sys
from shutil import which
from pathlib import Path
from importlib.metadata import version as get_pip_version, PackageNotFoundError
import yaml

# The environment file in the repository
ENVIRONMENT_FILE = Path(__file__).resolve().parents[2] / 'envs' / 'environment.yml'

# Comparison operators supported in version constraints
OPERATORS = {
    '>=': lambda a, b: a >= b,
    '<=': lambda a, b: a <= b,
    '==': lambda a, b: a == b,
    '!=': lambda a, b: a != b,
    '>': lambda a, b: a > b,
    '<': lambda a, b: a < b,
}


def version_key(version: str) -> list:
    """Split a version string into comparable pieces.
    Numbers are compared as numbers and letters (e.g. 'a', 'rc') go before numbers, so 1.0a1 < 1.0.
    """
    return [(1, int(piece)) if piece.isdigit() else (0, piece.lower())
        for piece in re.findall(r'\d+|[a-zA-Z]+', version)]


def compare_versions(a: str, b: str) -> int:
    """Return -1, 0 or 1 when version a is lower, equal or greater than version b."""
    key_a, key_b = version_key(a), version_key(b)
    # Pad the shorter version with zeros so 2.11 == 2.11.0
    size = max(len(key_a), len(key_b))
    key_a += [(1, 0)] * (size - len(key_a))
    key_b += [(1, 0)] * (size - len(key_b))
    return (key_a > key_b) - (key_a < key_b)


def parse_spec(spec: str) -> tuple[str, str]:
    """Split a dependency line into the package name and its version constraint.
    Channel prefixes (e.g. 'bioconda::') are removed.
    """
    spec = spec.split('::')[-1].strip()
    match = re.match(r'^([A-Za-z0-9_.\-]+)\s*(.*)$', spec)
    if not match: raise ValueError(f'Cannot read dependency "{spec}" in {ENVIRONMENT_FILE}')
    name, constraint = match.groups()
    return name.lower(), constraint.replace(' ', '')


def satisfies(installed_version: str, installed_build: str, constraint: str) -> bool:
    """Check if an installed version (and build) satisfies a conda version constraint."""
    if not constraint: return True
    if '|' in constraint: raise ValueError(f'OR constraints are not supported: "{constraint}"')
    # Conda style 'name=version' or 'name=version=build'
    if constraint.startswith('=') and not constraint.startswith('=='):
        pieces = constraint[1:].split('=')
        if len(pieces) > 2: raise ValueError(f'Cannot read constraint "{constraint}"')
        version = pieces[0].rstrip('*').rstrip('.')
        # '=2.7' accepts 2.7, 2.7.0, 2.7.1... but not 2.70
        if installed_version != version and not installed_version.startswith(version + '.'):
            return False
        if len(pieces) == 2 and pieces[1] != installed_build:
            return False
        return True
    # Constraints like '>=3.6.0,<4'
    for condition in constraint.split(','):
        match = re.match(r'^(>=|<=|==|!=|>|<)(.+)$', condition)
        if not match: raise ValueError(f'Cannot read constraint "{constraint}"')
        operator, version = match.groups()
        if not OPERATORS[operator](compare_versions(installed_version, version), 0):
            return False
    return True


def get_conda_packages() -> dict[str, tuple[str, str]]:
    """Get the installed conda packages with their version and build.
    Read the file names in conda-meta (name-version-build.json), which is much faster than calling conda.
    """
    conda_meta = Path(sys.prefix) / 'conda-meta'
    if not conda_meta.is_dir():
        raise SystemExit(f'The mwf command must run in a conda environment created from {ENVIRONMENT_FILE}\n'
            f'but the current python ({sys.executable}) is not in a conda environment')
    packages = {}
    for filename in conda_meta.glob('*.json'):
        name, version, build = filename.stem.rsplit('-', 2)
        packages[name.lower()] = (version, build)
    return packages


def check_dependencies():
    """Check every dependency in the environment file is installed with a valid version.
    Stop with a message listing all problems found otherwise.
    """
    # The environment file only exists when the workflow is installed from the repository (pip install -e .)
    # Packaged installs (e.g. bioconda or pip) do not ship it, but their dependencies were already resolved by the package manager
    if not ENVIRONMENT_FILE.exists(): return
    with open(ENVIRONMENT_FILE) as file:
        environment = yaml.safe_load(file)
    installed = get_conda_packages()
    problems = []
    for dependency in environment['dependencies']:
        # Pip dependencies come in a separated list
        if isinstance(dependency, dict):
            for spec in dependency.get('pip', []):
                name, constraint = parse_spec(spec)
                try:
                    version = get_pip_version(name)
                except PackageNotFoundError:
                    problems.append(f'  {name} (pip) is missing, required {constraint or "any version"}')
                    continue
                if not satisfies(version, '', constraint):
                    problems.append(f'  {name} (pip) is {version}, required {constraint}')
            continue
        name, constraint = parse_spec(dependency)
        if name not in installed:
            problems.append(f'  {name} is missing, required {constraint or "any version"}')
            continue
        version, build = installed[name]
        if not satisfies(version, build, constraint):
            problems.append(f'  {name} is {version} (build {build}), required {constraint}')
    if problems:
        # Suggest mamba if the user has it since it is much faster
        manager = 'mamba' if which('mamba') else 'conda'
        raise SystemExit(f'The environment in {sys.prefix} does not match {ENVIRONMENT_FILE}:\n'
            + '\n'.join(problems) + '\n'
            f'Update it with: {manager} env update --prefix {sys.prefix} --file {ENVIRONMENT_FILE}\n'
            f'If you are sure these versions are not a problem, skip this check by setting MWF_SKIPDEPS=1')
