#!/usr/bin/env python3
"""Measure a dataset's disk use, check its upload records, and list cleanup classes to delete.

Step 7 of the skill, clean up after the upload (see references/data-lifecycle.md). It never deletes
anything: it only prints sizes, checks and paths.

    dataset_inventory.py ROOT                           size per file name (digits -> N), largest first
    dataset_inventory.py ROOT --by ext                  the same by extension
    dataset_inventory.py ROOT --uploads                 check the upload records of every project
    dataset_inventory.py ROOT --classes classes.yaml    size per class, then the largest unclassified names
    dataset_inventory.py ROOT --classes classes.yaml --list CLASS [--relative] [-0] > CLASS.list
                                                        the paths of one class (or of `protected`)

A classes file (YAML, which needs PyYAML, or JSON) maps each class to glob patterns. A pattern without
'/' matches the file name at any depth; a pattern with '/' matches the path relative to ROOT, and its
'*' also crosses '/' ('IDRome_v4/*/traj.xtc' is any traj.xtc below IDRome_v4/). A pattern ending in
'/' is a folder: 'array_logs/' means everything below ROOT/array_logs. Patterns match files, never
dataset states: a class takes the files of every project, `done` or not.

    protected: ['run_logs/']            # optional: more files to keep, on top of PROTECTED below
    caches: ['IDRome_v4/*/mwf_*.out', '.*_offsets.npz']
    source_unused: ['IDRome.zip', 'emap.csv', '*.npy']

Files matching PROTECTED or the `protected` patterns (the dataset's bookkeeping) are never in a class:
a class that matches one is an error, and --list prints nothing until it is fixed. A file in several
classes counts in the first one and is reported, and so is a pattern that matches nothing. --list also
reports symlinks that would break: links outside the class that point to a listed file or to a folder
holding one. `--list protected` gives the bookkeeping files, for the backup before deleting anything.

--uploads finds the projects (folders with an inputs.yaml) and their upload records
(.loader_register.json, or .project_id from older loaders), and reports projects without one,
accessions shared by several projects, and records with no complete load (every load was run with a
flag that skips part of the data). It exits with 1 when it finds any of those.

On another host, without copying the script: ssh <host> python3 - <root> < dataset_inventory.py
"""
import argparse
import fnmatch
import json
import os
import re
import sys

# Bookkeeping that a cleanup must keep: upload records (the loader's register, and .project_id from
# older loaders), workflow caches and registers, inputs and metadata files, the dataset file and its
# README, the helper scripts and the dataset's own logs
PROTECTED = ['.loader_register.json', '.project_id', '.mwf_cache.json', '.register.json', 'inputs.yaml',
             'metadata.json', 'topology.json', '*.db', '*.sql', 'README.md',
             'utils/', 'run_logs/', 'array_logs/']
RESERVED = ('protected', 'unclassified')
# Loader flags that load only part of a project (see mddb_workflow/utils/loader.py)
PARTIAL_LOAD_FLAGS = {'--include', '--exclude', '--skip-trajectories', '--skip-files', '--skip-analyses',
                      '--md-directories'}


def walk(root: str):
    """Yield (path, is_symlink, size, bytes on disk) for every non-directory entry, without following links."""
    stack = [root]
    while stack:
        directory = stack.pop()
        try:
            entries = list(os.scandir(directory))
        except OSError as error:
            print(f'cannot read {directory}: {error}', file=sys.stderr)
            continue
        for entry in entries:
            try:
                if entry.is_dir(follow_symlinks=False):
                    stack.append(entry.path)
                    continue
                stat = entry.stat(follow_symlinks=False)
            except OSError as error:
                print(f'cannot read {entry.path}: {error}', file=sys.stderr)
                continue
            yield entry.path, entry.is_symlink(), stat.st_size, stat.st_blocks * 512


class Patterns:
    """Glob patterns matched against a file name (no '/') or the path relative to ROOT ('*' crosses '/')."""

    def __init__(self, patterns: list[str]):
        self.patterns = []
        for pattern in patterns:
            pattern = pattern.strip()
            if pattern.startswith('./'): pattern = pattern[2:]
            if pattern.endswith('/'): pattern += '*'
            self.patterns.append(pattern)
        self.regexes = {p: re.compile(fnmatch.translate(p)) for p in self.patterns}
        names = [fnmatch.translate(p) for p in self.patterns if '/' not in p]
        paths = [fnmatch.translate(p) for p in self.patterns if '/' in p]
        self.name_re = re.compile('|'.join(names)) if names else None
        self.path_re = re.compile('|'.join(paths)) if paths else None
        self.unused = set(self.patterns)

    def match(self, rel: str, name: str) -> bool:
        hit = bool((self.name_re and self.name_re.match(name)) or (self.path_re and self.path_re.match(rel)))
        # Remember which patterns matched something, to report the ones that never do
        if hit and self.unused:
            for pattern in list(self.unused):
                if self.regexes[pattern].match(rel if '/' in pattern else name):
                    self.unused.discard(pattern)
        return hit


def load_classes(path: str) -> tuple[dict, list]:
    """Read the classes file: ({class: [patterns]}, [extra protected patterns])."""
    with open(path) as file:
        if path.endswith(('.yaml', '.yml')):
            import yaml
            classes = yaml.safe_load(file)
        else:
            classes = json.load(file)
    if not isinstance(classes, dict):
        sys.exit(f'{path}: expected a mapping of class name -> list of glob patterns')
    if 'unclassified' in classes:
        sys.exit(f'{path}: "unclassified" is a reserved name')
    extra = classes.pop('protected', None) or []
    for name, patterns in [('protected', extra), *classes.items()]:
        if not isinstance(patterns, list) or not all(isinstance(p, str) and p.strip() for p in patterns):
            sys.exit(f'{path}: {name} must be a list of glob patterns (strings)')
        if name != 'protected' and not patterns:
            sys.exit(f'{path}: class {name} has no patterns')
    return classes, extra


def group_key(name: str, by: str) -> str:
    if by == 'ext':
        return os.path.splitext(name)[1] or '(no extension)'
    return re.sub(r'\d+', 'N', name)


def add(totals: dict, key: str, size: int, disk: int):
    row = totals.setdefault(key, [0, 0, 0])
    row[0] += size
    row[1] += disk
    row[2] += 1


def human(n: float) -> str:
    for unit, factor in (('TB', 1e12), ('GB', 1e9), ('MB', 1e6), ('KB', 1e3)):
        if n >= factor:
            return f'{n / factor:.1f} {unit}'
    return f'{n:.0f} B'


def plural(n: int, word: str) -> str:
    return f'{n:,} {word}' + ('' if n == 1 else 's')


def print_table(rows, label: str):
    print(f'{"size":>10} {"on disk":>10} {"count":>10}  {label}')
    for key, (size, disk, count) in rows:
        print(f'{human(size):>10} {human(disk):>10} {count:>10,}  {key}')


def check_uploads(projects: set, registers: dict, legacy: dict) -> bool:
    """Report projects without an upload record, shared accessions and records with no complete load."""
    missing = sorted(p for p in projects if p not in registers and p not in legacy)
    by_accession, partial, unreadable = {}, [], []
    for project, path in registers.items():
        try:
            with open(path) as file:
                entries = json.load(file)
            last = entries[-1]
            accession = last.get('accession') or last.get('projectId')
            complete = any(not PARTIAL_LOAD_FLAGS & set(str(e.get('command', '')).split()) for e in entries)
        except (OSError, ValueError, TypeError, AttributeError, IndexError) as error:
            unreadable.append(f'{project} ({error})')
            continue
        by_accession.setdefault(accession, []).append(project)
        if not complete: partial.append(project)
    for project, path in legacy.items():
        if project in registers: continue
        try:
            with open(path) as file:
                by_accession.setdefault(file.read().strip(), []).append(project)
        except OSError as error:
            unreadable.append(f'{project} ({error})')
    shared = {a: ps for a, ps in by_accession.items() if len(ps) > 1}
    outside = sorted((set(registers) | set(legacy)) - projects)
    legacy_only = sum(1 for p in legacy if p not in registers)
    print(f'\nUpload records: {plural(len(projects), "project")} (folders with an inputs.yaml), '
          f'{len(set(registers) | set(legacy)):,} with a record ({legacy_only:,} only .project_id), '
          f'{len(by_accession):,} different accessions')
    problems = [
        (missing, 'projects without an upload record'),
        ([f'{a}: {", ".join(ps)}' for a, ps in shared.items()], 'accessions in several projects'),
        (partial, 'records with no complete load (every load skipped part of the data)'),
        (unreadable, 'unreadable records'),
        (outside, 'records in folders without an inputs.yaml (not counted as projects)'),
    ]
    for i, (items, label) in enumerate(problems):
        if items:
            print(f'{"NOTE" if i == len(problems) - 1 else "PROBLEM"}: {label}: {len(items):,}, e.g.')
            for item in items[:5]: print(f'    {item}')
    return not (missing or shared or partial or unreadable)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('root', help='dataset root (or any folder holding a copy of it)')
    parser.add_argument('--by', choices=['name', 'ext'], default='name', help='group files by name or extension')
    parser.add_argument('--top', type=int, default=40, help='rows of the largest names (default 40)')
    parser.add_argument('--uploads', action='store_true', help='check the upload records of every project')
    parser.add_argument('--classes', help='YAML/JSON file mapping each class to glob patterns')
    parser.add_argument('--list', metavar='CLASS', help='print the paths of one class, or `protected` (needs --classes)')
    parser.add_argument('--relative', action='store_true', help='with --list, print paths relative to ROOT')
    parser.add_argument('-0', dest='null', action='store_true', help='with --list, separate paths with NUL')
    args = parser.parse_args()
    if args.list and not args.classes:
        parser.error('--list needs --classes')

    root = os.path.realpath(args.root)
    classes, extra_protected = load_classes(args.classes) if args.classes else ({}, [])
    if args.list and args.list not in [*classes, 'protected']:
        parser.error(f'--list {args.list}: not a class of {args.classes} ({", ".join(classes)}) nor protected')
    matchers = {name: Patterns(patterns) for name, patterns in classes.items()}
    protected = Patterns(PROTECTED + extra_protected)
    extra = Patterns(extra_protected)

    by_name, by_class = {}, {}
    conflicts, overlaps = {}, {}
    listed, links = [], []
    projects, registers, legacy = set(), {}, {}
    n_files = n_links = size_all = disk_all = 0
    for path, is_link, size, disk in walk(root):
        rel = path[len(root):].lstrip('/')
        name = os.path.basename(path)
        if is_link:
            n_links += 1
            if args.list: links.append(path)
        else:
            n_files += 1
            size_all += size
            disk_all += disk
        if args.uploads:
            folder = os.path.dirname(rel) or '.'
            if name == 'inputs.yaml': projects.add(folder)
            elif name == '.loader_register.json': registers[folder] = path
            elif name == '.project_id': legacy[folder] = path
        if not classes:
            if not is_link: add(by_name, group_key(name, args.by), size, disk)
            continue
        hits = [c for c, matcher in matchers.items() if matcher.match(rel, name)]
        extra.match(rel, name)
        if protected.match(rel, name):
            for c in hits: conflicts.setdefault(c, []).append(rel)
            group = 'protected'
        elif hits:
            if len(hits) > 1: overlaps.setdefault(' + '.join(hits), []).append(rel)
            group = hits[0]
        else:
            group = 'unclassified'
            if not is_link: add(by_name, group_key(name, args.by), size, disk)
        add(by_class, group, size, disk)
        if args.list and group == args.list:
            listed.append(path)

    unused = [(c, p) for c, matcher in [*matchers.items(), ('protected', extra)] for p in sorted(matcher.unused)]

    if args.list:
        if conflicts:
            for c, rels in conflicts.items():
                print(f'ERROR: class {c} matches {plural(len(rels), "protected file")}, e.g. {rels[0]}', file=sys.stderr)
            sys.exit('Nothing listed: remove those patterns from the classes file first.')
        sep = '\0' if args.null else '\n'
        start = len(root) + 1 if args.relative else 0
        sys.stdout.write(''.join(p[start:] + sep for p in listed))
        size, disk, _ = by_class.get(args.list, (0, 0, 0))
        print(f'{args.list}: {plural(len(listed), "path")}, {human(size)} ({human(disk)} on disk)', file=sys.stderr)
        for c, pattern in unused:
            if c == args.list: print(f'WARNING: pattern {pattern!r} of {c} matches nothing', file=sys.stderr)
        # Links outside the class that point to a listed file, or to a folder holding one
        listed_set = set(listed)
        folders = set()
        for parent in {os.path.dirname(p) for p in listed}:
            while parent.startswith(root) and parent not in folders:
                folders.add(parent)
                parent = os.path.dirname(parent)
        targets = listed_set | folders
        broken = [l for l in links if l not in listed_set and os.path.realpath(l) in targets]
        if broken:
            print(f'WARNING: {plural(len(broken), "symlink")} outside {args.list} would break (they point to a '
                  f'listed file or a folder holding one), e.g. {broken[0]}', file=sys.stderr)
        return

    print(f'{root}: {n_files:,} files, {human(size_all)} ({human(disk_all)} on disk), and {n_links:,} symlinks\n')
    if classes:
        order = [c for c in [*classes, *RESERVED] if c in by_class]
        print_table([(c, by_class[c]) for c in order], 'class')
        for c, rels in conflicts.items():
            print(f'\nERROR: class {c} matches {plural(len(rels), "protected file")}, e.g. {rels[0]}')
        for pair, rels in overlaps.items():
            print(f'\nWARNING: {plural(len(rels), "file")} in {pair} (counted in the first), e.g. {rels[0]}')
        for c, pattern in unused:
            print(f'\nWARNING: pattern {pattern!r} of {c} matches nothing')
        if by_name:
            print('\nLargest unclassified names:')
    if by_name:
        rows = sorted(by_name.items(), key=lambda kv: -kv[1][1])[:args.top]
        print_table(rows, 'extension' if args.by == 'ext' else 'name (digits -> N)')
    if args.uploads and not check_uploads(projects, registers, legacy):
        sys.exit(1)


if __name__ == '__main__':
    main()
