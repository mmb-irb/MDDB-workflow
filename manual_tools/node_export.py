"""Export a project already uploaded to an MDDB node into a directory the loader can load as is.

PROTOTYPE (investigation of node-to-node migration). Only reads from the node REST API.

The export is a byte-for-byte copy of what the node stores, laid out so the existing loader `load` command
uploads it without re-deriving anything:
  - project `metadata.json`, `topology.json` and reference files (pdb_references.json, collections.json...);
  - one directory per MD with the MD document (minus its file and analysis lists) as `metadata.json`;
  - every analysis as `mda.<name>.json`, holding the stored value;
  - every GridFS file as `mdf.<name>` (the loader strips the prefix), with `mdf.<name>.meta.json` carrying the
    file document metadata (e.g. frames and atoms of .bin files) besides project and md.
    `trajectory.bin` is copied as stored, so no `gmx dump` runs again.
Exceptions to the `mdf.` prefix: `structure.*` and `topology.*` keep their name (the loader needs them as such).
The loader requires a `trajectory.xtc` in every MD directory of a new project, so a hard link to
`mdf.trajectory.xtc` is left there: load with `--skip-trajectories` so it is not parsed or uploaded twice.

No inputs file is written, and the source accession is not carried in any loaded file, so the target node issues
a new accession. Top-level project fields the loader does not take from files (published, mdref, dates) are left
in `.export_source.json` together with the source file sizes, to check the copy afterwards.

Usage:
    python -I node_export.py <api_url> <accession> <output_dir> [--metadata-only]
e.g. python -I node_export.py https://irb-dev.mddbr.eu/api/ A05DO exports/A05DO
"""
import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

CHUNK_SIZE = 4 * 1024 * 1024
RETRIES = 3
# Map API reference types to the loader reference file names
REFERENCE_FILES = {
    'proteins': 'protein_references.json',
    'inchikeys': 'inchikey_references.json',
    'pdbs': 'pdb_references.json',
    'chains': 'chains.json',
    'collections': 'collections.json',
}
# File names the loader must find without the 'mdf.' prefix
KEEP_NAME_PREFIXES = ('structure.', 'topology.')
# Project document fields which are rebuilt by the loader and must not go to the MD metadata files
MD_LIST_FIELDS = ('files', 'analyses')


def request(url: str, output_path: str = None, expected_size: int = None):
    """GET a URL, retried. Return the parsed JSON, or stream it to a file and check its size."""
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                if output_path is None:
                    return json.loads(response.read())
                written = 0
                with open(output_path + '.part', 'wb') as file:
                    while chunk := response.read(CHUNK_SIZE):
                        file.write(chunk)
                        written += len(chunk)
                if expected_size is not None and written != expected_size:
                    raise IOError(f'{url}: got {written} B, expected {expected_size} B')
                os.replace(output_path + '.part', output_path)
                return written
        except (urllib.error.URLError, IOError, TimeoutError) as error:
            if attempt == RETRIES:
                raise
            print(f'  retry {attempt} after error: {error}', file=sys.stderr)
            time.sleep(2 ** attempt)


def save_json(data, path: str):
    with open(path, 'w') as file:
        json.dump(data, file)


def local_filename(name: str) -> str:
    return name if name.startswith(KEEP_NAME_PREFIXES) else f'mdf.{name}'


def export_project(api_url: str, accession: str, output_dir: str, metadata_only: bool = False):
    base = api_url.rstrip('/') + '/rest/v1/projects/'
    quote = lambda value: urllib.parse.quote(str(value), safe='')
    project = request(f'{base}{quote(accession)}?raw')
    descriptors = request(f'{base}{quote(accession)}/filenotes')
    topology = request(f'{base}{quote(accession)}/topology')
    references = request(f'{base}{quote(accession)}/references')
    os.makedirs(output_dir, exist_ok=False)

    # Project level documents
    save_json(project['metadata'], os.path.join(output_dir, 'metadata.json'))
    for field in ('_id', 'project'):
        topology.pop(field, None)
    save_json(topology, os.path.join(output_dir, 'topology.json'))
    by_type = {}
    for reference in references:
        reference = dict(reference)
        by_type.setdefault(reference.pop('ref_type'), []).append(reference)
    for ref_type, docs in by_type.items():
        if ref_type not in REFERENCE_FILES:
            raise ValueError(f'Unknown reference type {ref_type}')
        save_json(docs, os.path.join(output_dir, REFERENCE_FILES[ref_type]))

    # MD directories, named after their index so names with odd characters do not matter
    # The loader finds MD directories by their '.register.json' file
    md_dirs = {}
    for md_index, md in enumerate(project['mds']):
        if md.get('removed'):
            continue
        md_dir = os.path.join(output_dir, f'replica_{md_index + 1}')
        os.makedirs(md_dir)
        save_json({'exported_from': f'{api_url} {accession}.{md_index + 1}'}, os.path.join(md_dir, '.register.json'))
        save_json({k: v for k, v in md.items() if k not in MD_LIST_FIELDS}, os.path.join(md_dir, 'metadata.json'))
        md_dirs[md_index] = md_dir

    # Analyses: the API serves the stored value, which is what the loader expects in mda.<name>.json
    # Project level analyses (no MD) are asked to the project and MD analyses to '<accession>.<md number>'
    analyses_to_get = [(None, a['name']) for a in project.get('analyses', [])]
    for md_index in md_dirs:
        analyses_to_get += [(md_index, a['name']) for a in project['mds'][md_index].get('analyses', [])]
    for md_index, name in analyses_to_get:
        target = accession if md_index is None else f'{accession}.{md_index + 1}'
        value = request(f'{base}{quote(target)}/analyses/{quote(name)}')
        folder = output_dir if md_index is None else md_dirs[md_index]
        save_json(value, os.path.join(folder, f'mda.{name}.json'))

    # Files, as listed by their GridFS documents, which is what really exists
    # Names listed in the project document without a GridFS document are reported, not copied
    listed = {(None, f['name']) for f in project.get('files', [])}
    for md_index in md_dirs:
        listed |= {(md_index, f['name']) for f in project['mds'][md_index].get('files', [])}
    stored = set()
    sizes = {}
    for descriptor in descriptors:
        metadata = descriptor.get('metadata', {})
        md_index = metadata.get('md')
        name = descriptor['filename']
        stored.add((md_index, name))
        folder = output_dir if md_index is None else md_dirs.get(md_index)
        if folder is None:
            print(f'  skipping {name}: it belongs to a removed or unknown MD {md_index}', file=sys.stderr)
            continue
        filename = local_filename(name)
        extra_metadata = {k: v for k, v in metadata.items() if k not in ('project', 'md')}
        if extra_metadata:
            save_json(extra_metadata, os.path.join(folder, filename + '.meta.json'))
        sizes[os.path.relpath(os.path.join(folder, filename), output_dir)] = descriptor['length']
        if metadata_only:
            continue
        target = accession if md_index is None else f'{accession}.{md_index + 1}'
        print(f'  {target} {name} ({descriptor["length"] / 1e6:.1f} MB)')
        request(f'{base}{quote(target)}/files/{quote(name)}', os.path.join(folder, filename), descriptor['length'])
        if name == 'trajectory.xtc' and md_index is not None:
            os.link(os.path.join(folder, filename), os.path.join(folder, 'trajectory.xtc'))
    dangling = sorted(listed - stored, key=str)

    save_json({
        'api_url': api_url,
        'accession': accession,
        'source_id': project['_id'],
        'not_loaded_fields': {k: project.get(k) for k in ('published', 'mdref', 'creationDate', 'updateDate')},
        'file_sizes': sizes,
        'dangling_file_entries': dangling,
        'analyses': [list(a) for a in analyses_to_get],
    }, os.path.join(output_dir, '.export_source.json'))
    total = sum(sizes.values())
    print(f'{accession}: {len(sizes)} files ({total / 1e9:.2f} GB), {len(analyses_to_get)} analyses, '
          f'{len(by_type)} reference types, {len(md_dirs)} MDs, dangling file entries: {dangling}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('api_url')
    parser.add_argument('accession')
    parser.add_argument('output_dir')
    parser.add_argument('--metadata-only', action='store_true', help='write everything but the file contents')
    args = parser.parse_args()
    export_project(args.api_url, args.accession, args.output_dir, args.metadata_only)


if __name__ == '__main__':
    main()
