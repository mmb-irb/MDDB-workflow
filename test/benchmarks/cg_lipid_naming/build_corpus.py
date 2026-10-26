"""Build the labelled corpus of coarse grained residues used to benchmark the recognition of lipids by bead names.

The parameter sets of Martini and other coarse grained force fields are downloaded into a cache directory (git
repositories at pinned commits, cgmartini and SIRAH files), their topologies are parsed as text and every residue is
labelled as a fatty lipid ('fatty'), a sterol ('sterol') or not a lipid (None) from the file or section it comes from.
Downloaded files are only read as text, archives are read member by member without extracting them.

The unique (residue name, bead names, label) cases are written to corpus.json, with the bead names separated by spaces,
the class of molecule, the source groups where they are found and the first file where they are found.
See README.md for the results and test_cg_lipid_naming.py for the evaluation of the rule in mddb_workflow.

Usage:
    python test/benchmarks/cg_lipid_naming/build_corpus.py [--cache DIR] [--output corpus.json]
"""
import argparse
import glob
import json
import os
import re
import subprocess
import tarfile
import zipfile
import requests
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
# Git repositories of parameter sets, at the commits of the benchmark
REPOSITORIES = {
    'M3-Lipid-Parameters': ('Martini-Force-Field-Initiative/M3-Lipid-Parameters', '2121527045696345069deb4d3df754724888a4f9'),
    'M3-Sterol-Parameters': ('Martini-Force-Field-Initiative/M3-Sterol-Parameters', '421efdaecfe54be5efd2076f9c7c706e79746ecb'),
    'martini-forcefields': ('marrink-lab/martini-forcefields', '784591ebdc91d762ed4df986c4650546c938f776'),
    'M3_Ionizable_Lipids': ('Martini-Force-Field-Initiative/M3_Ionizable_Lipids', '9954b0f85318045d883f9d12db4bff3749c5d7d5'),
    'Martini3-small-molecules': ('ricalessandri/Martini3-small-molecules', '9a630c18934be270b4f77189722650f4397e4ef2'),
    'polyply_1.0': ('marrink-lab/polyply_1.0', '60cc6aa3ad057126e586f4cb8919b2e0ac069ab1'),
    'vermouth-martinize': ('marrink-lab/vermouth-martinize', 'a9b62eed6eec028c486ab6369a13cdd9c31f4561'),
    'Insane': ('Tsjerk/Insane', 'e23a501036c683f2652fe07a66ba941195d7b8f6'),
    'TS2CG1.1': ('marrink-lab/TS2CG1.1', 'a5ee11888e4cda110bc69410df7070ef66090adc'),
    'spica-tools': ('SPICA-group/spica-tools', '285a1ffa900ce5b6b7fa5b95277300031572943b'),
}
# Martini 2 lipidome page, which links every lipid ITP
M2_LIPIDOME_PAGE = 'https://cgmartini.nl/docs/downloads/force-field-parameters/martini2/lipidome.html'
# Other Martini 2 parameter files (non-lipids and a few lipids), relative to the cgmartini library
M2_LIBRARY = 'https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini2'
M2_OTHER_FILES = [
    'carbohydrates/martini_v2.0_sugars.itp', 'others/martini_v2.0_other.itp', 'others/martini_v2.0_phytochemicals.itp',
    'others/nucleotides.itp', 'others/martini_v2.0_CHOL_02-optLINCS.itp', 'polymers/martini_v2.0_polymers.itp',
    'polymers/peg/martini_v2.0_PEO_PS_CNP.itp', 'polymers/peg/martini_v2.2refP_PEO.itp', 'polymers/pss/martini_v2.2P_PS.itp',
    'solvents/martini_v2.0_solvents.itp', 'aminoacids/martini_v2.2_aminoacids.itp', 'aminoacids/martini_v2.2P_aminoacids.itp',
    'others/topologies.zip', 'others/BCR_correct.zip', 'dry_martini/Dry_Martini_v2.1_20160523.zip',
    'nucleic_acids/dna/martini-dna-150909.tar', 'nucleic_acids/rna/na-tutorials_20170815.tar',
]
# SIRAH force field for GROMACS
SIRAH_ARCHIVE = 'https://github.com/SIRAHFF/documentation/releases/download/GROMACS/sirah_x2.2_20-07.ff.tar.gz'
# Residue names of SPICA lipids and sterols (other SPICA residues are not lipids)
SPICA_LIPIDS = {'CER1', 'CER2', 'DAPC', 'DLPC', 'DMPC', 'DOPC', 'POPE', 'DOPE', 'DPPC', 'DSPC', 'DLIPC', 'DLIP', 'POPC',
    'DOPG', 'DMPG', 'POPG', 'SPM2', 'SSM', 'PSM', 'ASM', 'LSM', 'BSM', 'POPI', 'POPS', 'DOPS', 'BMGP'}
SPICA_STEROLS = {'CHL1', 'CLR', 'CLOL'}
# Topology sections which are not residues in SIRAH residue topology files
SIRAH_SECTIONS = {'atoms', 'bonds', 'angles', 'dihedrals', 'impropers', 'exclusions', 'cmap', 'bondedtypes'}


def download(url: str, path: str) -> str:
    """Download a file unless it is already in the cache."""
    if not os.path.exists(path):
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'wb') as handle:
            handle.write(response.content)
    return path


def checkout(cache: str, name: str) -> str:
    """Get a git repository at its pinned commit, without history."""
    repository, commit = REPOSITORIES[name]
    path = os.path.join(cache, name)
    if not os.path.exists(os.path.join(path, '.git')):
        subprocess.run(['git', 'init', '-q', path], check=True)
        subprocess.run(['git', '-C', path, 'fetch', '-q', '--depth', '1', f'https://github.com/{repository}.git', commit], check=True)
        subprocess.run(['git', '-C', path, 'checkout', '-q', 'FETCH_HEAD'], check=True)
    return path


def parse_itp(lines: list[str]) -> list[dict]:
    """Get the molecules of a GROMACS topology (or a vermouth/polyply .ff file) with their residues and bead names."""
    molecules = []
    section = None
    current = None
    for raw in lines:
        line = raw.split(';')[0].strip()
        if not line or line.startswith('#'): continue
        match = re.match(r'^\[\s*(\w+)\s*\]', line)
        if match:
            section = match.group(1).lower()
            continue
        if section == 'moleculetype':
            current = {'name': line.split()[0], 'residues': OrderedDict()}
            molecules.append(current)
            section = None
        elif section == 'atoms' and current is not None:
            fields = line.split()
            if len(fields) < 5 or not fields[2].lstrip('-').isdigit(): continue
            current['residues'].setdefault((int(fields[2]), fields[3]), []).append(fields[4])
    return molecules


def read_lines(path: str) -> list[str]:
    with open(path, errors='replace') as handle:
        return handle.read().splitlines()


def archive_members(path: str, pattern: str) -> list[tuple[str, list[str]]]:
    """Get the text lines of the archive members whose name matches a pattern, without extracting them."""
    members = []
    if path.endswith('.zip'):
        with zipfile.ZipFile(path) as archive:
            for name in sorted(archive.namelist()):
                if re.search(pattern, name):
                    members.append((name, archive.read(name).decode(errors='replace').splitlines()))
    else:
        with tarfile.open(path) as archive:
            for member in sorted(archive.getmembers(), key=lambda member: member.name):
                if member.isfile() and re.search(pattern, member.name):
                    members.append((member.name, archive.extractfile(member).read().decode(errors='replace').splitlines()))
    return members


class Corpus:
    """Labelled residues by (residue name, bead names, label), with the groups and the first file where they are found."""

    def __init__(self):
        self.cases = OrderedDict()

    def add(self, group: str, source: str, resname: str, beads: list[str], expected: str | None, label: str, note: str = ''):
        key = (resname, tuple(beads), expected)
        if key not in self.cases:
            self.cases[key] = {'resname': resname, 'beads': ' '.join(beads), 'expected': expected, 'class': label,
                'groups': [], 'source': source}
            if note: self.cases[key]['note'] = note
        if group not in self.cases[key]['groups']: self.cases[key]['groups'].append(group)

    def add_itp(self, group: str, source: str, lines: list[str], classify: callable):
        """Add the residues of a topology, labelled by a function of the residue name, its bead names and whether the
        molecule has several residues. Residues labelled with no class are skipped.
        """
        for molecule in parse_itp(lines):
            multi = len(molecule['residues']) > 1
            for (_, resname), beads in molecule['residues'].items():
                expected, label, note = classify(resname, beads, multi)
                if label is None: continue
                self.add(group, source, resname, beads, expected, label, note)


def m3_lipidome(filename: str) -> callable:
    def classify(resname, beads, multi):
        if 'ions' in filename or 'solvents' in filename: return None, 'solvent/ion', ''
        if 'ffbonded' in filename or filename == 'martini_v3.0.0.itp': return None, None, ''
        if 'hydrocarbons' in filename: return None, 'alkane (no head)', 'pure hydrocarbon oil, no polar bead'
        if 'sterols' in filename: return 'sterol', 'sterol', ''
        if resname.startswith('INO'): return None, 'free inositol phosphate', ''
        return 'fatty', re.sub(r'martini_v3\.0\.0_|_v\d+\.itp', '', filename), ''
    return classify


def m3_official(filename: str) -> callable:
    def classify(resname, beads, multi):
        if 'phospholipids' in filename: return 'fatty', 'phospholipid (M3 v1)', ''
        if filename == 'martini_v3.0.0.itp': return None, None, ''
        return None, re.sub(r'martini_v3\.0\.0_|_v\d+\.itp|\.itp', '', filename), ''
    return classify


def m3_ionizable(relative: str) -> callable:
    def classify(resname, beads, multi):
        # Atomistic mapping topologies
        if '/Mapping/' in relative: return None, None, ''
        if 'Fragments' in relative: return None, 'IL fragment', 'parameterisation fragment, not a lipid'
        if resname in ('PEO', 'OHter'): return None, 'PEG monomer of PEG-lipid', 'separate residue of a PEG-lipid molecule'
        if 'Sterols' in relative or 'sterols' in os.path.basename(relative): return 'sterol', 'sterol', ''
        if 'Literature_known_lipids' in relative or 'Ionizable_lipid_library' in relative:
            if 'ffbonded' in relative: return None, None, ''
            return 'fatty', 'ionizable/PEG lipid', ''
        return None, None, ''
    return classify


def m2_lipidome(filename: str) -> callable:
    category = filename.split('__')[0]
    classes = {'pa': 'PA', 'pc': 'PC', 'pe': 'PE', 'pg': 'PG', 'ps': 'PS', 'pi': 'PI/PIP', 'sm': 'SM', 'cer': 'ceramide',
        'gl': 'glycoglycerolipid', 'gs': 'glycosphingolipid', 'glycerols': 'DAG/TAG'}

    def classify(resname, beads, multi):
        if category == 'sterols':
            if resname in ('HOPR', 'HHOP', 'HDPT', 'HBHT'): return 'sterol', 'hopanoid', ''
            if resname == 'CHYO': return 'sterol', 'sterol ester', 'cholesteryl oleate, core lipid'
            return 'sterol', 'sterol', ''
        if category == 'surfactants': return 'fatty', 'detergent', 'single-tail detergent'
        if category == 'fa': return 'fatty', 'fatty acid (M2)', ''
        if category == 'lpc': return 'fatty', 'lyso-PC', ''
        if category == 'others':
            if resname.startswith('CDL'): return 'fatty', 'cardiolipin', ''
            if resname.startswith('BOL'): return 'fatty', 'bolalipid', ''
            if 'TAP' in resname: return 'fatty', 'cationic lipid', ''
            if multi: return None, 'lipid A / LPS residue', 'one residue of a multi-residue lipid A molecule'
            return 'fatty', 'LPS', 'single-residue LPS'
        return 'fatty', classes.get(category, category), ''
    return classify


def m2_other(filename: str) -> callable:
    def classify(resname, beads, multi):
        if 'CHOL' in filename or 'cholesterol' in filename: return 'sterol', 'sterol', ''
        if 'dry_martini_v2.1_lipids' in filename: return 'fatty', 'dry Martini lipid', ''
        if filename.endswith('martini_v2.0_polymers.itp') and resname == 'PEL': return 'fatty', 'PEG-lipid (M2)', ''
        if 'solvents' in filename and resname == 'TO': return 'fatty', 'TAG (M2 solvents file)', ''
        if 'other' in filename and resname.startswith(('DYE', 'ALT')):
            return None, 'membrane dye', 'voltage-sensitive styryl dye (di-n-ASP), not a lipid'
        if 'heme' in filename: return None, 'heme', ''
        if any(key in filename for key in ('ubiquinone', 'plastoquinone', 'betacarot', 'chlorophyll', 'pheophytin')):
            return None, 'cofactor/pigment', ''
        if 'dna' in filename.lower() or 'Nucleic' in filename: return None, 'DNA (M2)', ''
        if 'aminoacids' in filename: return None, None, ''
        if filename.endswith('_PS.itp'): return None, 'polystyrene + polarizable water (M2 2.2P)', ''
        return None, re.sub(r'.*martini_v2\.\w*?_|\.itp', '', filename), ''
    return classify


def polyply_vermouth(path: str) -> callable:
    def classify(resname, beads, multi):
        if 'aminoacids' in path:
            if resname in ('GLYM', 'CYSP', 'CYSF', 'CYSG'): return None, 'lipidated amino acid', 'protein residue with a lipid anchor'
            return None, 'protein (Martini)', ''
        return None, 'polymer/other (polyply/vermouth)', ''
    return classify


def build_corpus(cache: str) -> Corpus:
    corpus = Corpus()
    paths = {name: checkout(cache, name) for name in REPOSITORIES}

    def add_files(group, root, patterns, classify, prefix):
        files = sorted(set(path for pattern in patterns for path in glob.glob(os.path.join(root, pattern), recursive=True)))
        for path in files:
            relative = os.path.relpath(path, root)
            corpus.add_itp(group, f'{prefix}:{relative}', read_lines(path), classify(relative))

    # Martini 3
    lipidome = paths['M3-Lipid-Parameters']
    add_files('M3 lipidome', lipidome, ['ITPs/*.itp'], lambda relative: m3_lipidome(os.path.basename(relative)), 'M3-Lipid-Parameters')
    add_files('M3 sterols', paths['M3-Sterol-Parameters'], ['martini_v3.0_sterols_v1.0.itp'],
        lambda relative: lambda *_: ('sterol', 'sterol', ''), 'M3-Sterol-Parameters')
    official = paths['martini-forcefields']
    for path in sorted(glob.glob(f'{official}/martini_forcefields/regular/v3.0.0/gmx_files*/*.itp')):
        filename = os.path.basename(path)
        group = 'M3 v1 phospholipids' if 'phospholipids' in filename else 'M3 official non-lipids'
        corpus.add_itp(group, 'martini-forcefields:' + os.path.relpath(path, official), read_lines(path), m3_official(filename))
    # Coarse grained models only, not the all atom (OPLS) models of the library
    add_files('M3 small-molecule library', paths['Martini3-small-molecules'], ['models/*.itp', 'models/itps/**/*.itp'],
        lambda relative: lambda *_: (None, 'small molecule', ''), 'Martini3-small-molecules')
    ionizable = paths['M3_Ionizable_Lipids']
    add_files('M3 ionizable-lipid repo', ionizable, ['Lipid_parameters/**/*.itp'], m3_ionizable, 'M3_Ionizable_Lipids')
    add_files('M3 ionizable-lipid repo', ionizable,
        ['Case_studies/Lipid_Nanoparticle_construction/Bleb_structure/Building_notebooks/ITP/martini_v3.0_sterols_v1.0.itp'],
        lambda relative: lambda *_: ('sterol', 'sterol', ''), 'M3_Ionizable_Lipids')
    add_files('M3 ionizable-lipid repo', ionizable, ['support/RNA.itp'],
        lambda relative: lambda *_: (None, 'RNA (M3)', ''), 'M3_Ionizable_Lipids')
    ts2cg = paths['TS2CG1.1']
    add_files('TS2CG', ts2cg, ['Tutorials/files/itp/martini3/martini_v3.0_CDLs.itp'],
        lambda relative: lambda *_: ('fatty', 'cardiolipin (M3 dev)', ''), 'TS2CG1.1')

    # Martini 2 lipidome, every lipid ITP linked from its page
    page = requests.get(M2_LIPIDOME_PAGE, timeout=120)
    page.raise_for_status()
    urls = sorted(set(re.findall(r'href="(https://cgmartini-library[^"]+\.itp)"', page.text)))
    for url in urls:
        filename = url.split('/lipidome/')[-1].replace('/', '__')
        path = download(url, os.path.join(cache, 'm2_lipidome', filename))
        corpus.add_itp('M2 lipidome', 'cgmartini M2 lipidome:' + filename, read_lines(path), m2_lipidome(filename))

    # Other Martini 2 files and archives
    for file in M2_OTHER_FILES:
        path = download(f'{M2_LIBRARY}/{file}', os.path.join(cache, 'm2_other', file.replace('/', '_')))
        source = 'cgmartini M2:' + os.path.basename(path)
        if path.endswith('.itp'):
            corpus.add_itp('M2 other', source, read_lines(path), m2_other(os.path.basename(path)))
            continue
        for member, lines in archive_members(path, r'\.itp$'):
            if 'MACOSX' in member or 'gromos' in member or 'tip3p' in member: continue
            corpus.add_itp('M2 other', f'{source} :: {member}', lines, m2_other(os.path.basename(member)))

    # polyply and vermouth libraries
    add_files('polyply/vermouth', paths['polyply_1.0'], ['polyply/data/martini[23]/*.ff', 'polyply/data/martini[23]/*.itp'],
        polyply_vermouth, 'polyply_1.0')
    add_files('polyply/vermouth', paths['vermouth-martinize'],
        ['vermouth/data/force_fields/martini3001/*.ff', 'vermouth/data/force_fields/martini22/*.ff'], polyply_vermouth, 'vermouth-martinize')

    # insane lipid templates
    block = None
    for line in read_lines(os.path.join(paths['Insane'], 'insane', 'lipids.dat')):
        line = line.split(';')[0].rstrip()
        if not line.strip(): continue
        match = re.match(r'^\[\s*(.+?)\s*\]', line)
        if match:
            block = match.group(1)
            continue
        fields = line.split()
        if not re.match(r'^[A-Za-z0-9]+\.[A-Za-z][A-Za-z0-9_]*$', fields[0]) or len(fields) <= 2: continue
        beads = [field for field in fields[2:] if field != '-']
        if block in ('quinones', 'hexaperihexabenzocoronene'): expected, label = None, block
        elif block == 'sterols': expected, label = 'sterol', 'sterol'
        elif block == 'hopanoids': expected, label = 'sterol', 'hopanoid'
        elif block == 'retinylesters': expected, label = None, 'retinyl ester'
        else: expected, label = 'fatty', block
        corpus.add('insane', f'insane lipids.dat [{block}]', fields[0].split('.')[1], beads, expected, label)

    # TS2CG Martini 3 lipid library
    templates, current = {}, None
    for line in read_lines(os.path.join(ts2cg, 'Tutorials', 'files', 'Martini3.LIB')):
        line = line.split(';')[0].strip()
        match = re.match(r'^\[(.+)\]', line)
        if match:
            current = match.group(1).strip()
            templates[current] = []
            continue
        fields = line.split()
        if current and len(fields) >= 2 and fields[0].isdigit(): templates[current].append(fields[1])
    for name, beads in templates.items():
        corpus.add('TS2CG', 'TS2CG1.1 Martini3.LIB', name, beads, 'sterol' if name == 'CHOL' else 'fatty', 'TS2CG template')

    # SPICA
    with open(os.path.join(paths['spica-tools'], 'src', 'spica_top.json')) as handle:
        spica = json.load(handle)['topo']
    for name, topology in spica.items():
        expected = 'fatty' if name in SPICA_LIPIDS else 'sterol' if name in SPICA_STEROLS else None
        corpus.add('SPICA', 'SPICA spica_top.json', name, topology['name'], expected, 'SPICA ' + ('lipid' if expected else 'other'))

    # SIRAH residue topologies of lipids and the S-palmitoyl cysteine
    archive = download(SIRAH_ARCHIVE, os.path.join(cache, 'sirah', os.path.basename(SIRAH_ARCHIVE)))
    for member, lines in archive_members(archive, r'/(lipid[^/]*|aminoacids)\.rtp$'):
        filename = os.path.basename(member)
        residues, residue, section = OrderedDict(), None, None
        for line in lines:
            line = line.split(';')[0].strip()
            if not line: continue
            match = re.match(r'^\[\s*(\S+)\s*\]', line)
            if match:
                tag = match.group(1)
                if tag in SIRAH_SECTIONS: section = tag if tag != 'bondedtypes' else None
                else: residue, section = tag, None; residues[residue] = []
                continue
            if section == 'atoms' and residue: residues[residue].append(line.split()[0])
        for resname, beads in residues.items():
            if filename == 'aminoacids.rtp':
                if resname == 'sCp': corpus.add('SIRAH', 'SIRAH aminoacids.rtp', resname, beads, None, 'lipidated amino acid (SIRAH)')
                continue
            fragment = resname.startswith('x')
            corpus.add('SIRAH', 'SIRAH ' + filename, resname, beads, None if fragment else 'fatty',
                'SIRAH lipid fragment' if fragment else 'SIRAH lipid')
    return corpus


def write_corpus(corpus: Corpus, output: str):
    """Write the cases one per line, so changes are easy to review."""
    lines = [json.dumps(case) for case in corpus.cases.values()]
    with open(output, 'w') as handle:
        handle.write('[\n' + ',\n'.join(lines) + '\n]\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--cache', default=os.path.join(HERE, 'downloads'), help='Directory for the downloaded parameter sets')
    parser.add_argument('--output', default=os.path.join(HERE, 'corpus.json'), help='Output corpus file')
    args = parser.parse_args()
    corpus = build_corpus(args.cache)
    write_corpus(corpus, args.output)
    expected = [case['expected'] for case in corpus.cases.values()]
    print(f'{len(expected)} cases: {sum(1 for label in expected if label)} lipids and {expected.count(None)} non-lipids')
