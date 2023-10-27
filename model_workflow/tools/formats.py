# Get a file extension
def get_file_extension (filename : str) -> str:
    extension = filename.split('.')[-1]
    if extension == filename:
        raise Exception('File "' + filename + '" has not extension and thus format cannot be guessed')
    return extension

# Get a file format
# Note that some formats have different filename tails
standard_formats = {
    # Topologies
    'top': 'top',
    'psf': 'psf',
    'prmtop': 'prmtop',
    'prm7': 'prmtop',
    # Structures
    'pdb': 'pdb',
    'gro': 'gro',
    # Trajectories
    'xtc': 'xtc',
    'trr': 'trr',
    'dcd': 'dcd',
    'nc': 'nc',
    'cdf': 'nc',
    'netcdf': 'nc',
    'crd': 'crd',
    'mdcrd': 'crd',
    'trj': 'crd'
}
def get_file_standard_format (filename : str) -> str:
    extension = get_file_extension(filename)
    standard_format = standard_formats.get(extension, None)
    if not standard_format:
        raise Exception('Not recognized format extension "' + extension + '" from file "' + filename + '"')
    return standard_format

# Topology file formats

def is_pdb (filename : str) -> bool:
    return filename[-4:] == '.pdb'

def is_psf (filename : str) -> bool:
    return filename[-4:] == '.psf'

def is_tpr (filename : str) -> bool:
    return filename[-4:] == '.tpr'

def is_gro (filename : str) -> bool:
    return filename[-4:] == '.gro'

# Trajectory file formats

def is_xtc (filename : str) -> bool:
    return filename[-4:] == '.xtc'

def is_dcd (filename : str) -> bool:
    return filename[-4:] == '.dcd'

def is_netcdf (filename : str) -> bool:
    return filename[-3:] == '.nc'

def are_xtc (filenames : list) -> bool:
    return all([ is_xtc(filename) for filename in filenames ])

def are_dcd (filenames : list) -> bool:
    return all([ is_dcd(filename) for filename in filenames ])

def are_netcdf (filenames : list) -> bool:
    return all([ is_netcdf(filename) for filename in filenames ])

# Charges file formats
raw_charges_filename = 'charges.txt'
def is_raw (filename : str) -> bool:
    return filename == raw_charges_filename

def is_prmtop (filename : str) -> bool:
    return filename[-7:] == '.prmtop'

def is_top (filename : str) -> bool:
    return filename[-4:] == '.top'

def is_psf (filename : str) -> bool:
    return filename[-4:] == '.psf'

# Extra formats logic

# Check if a file may be read by pytraj according to its format
def is_pytraj_supported (filename : str) -> bool:
    return is_prmtop(filename) or is_top(filename) or is_psf(filename)

# From GitHub:
# ParmFormatDict = {
#     "AMBERPARM": AMBERPARM,
#     "PDBFILE": PDBFILEPARM,
#     "MOL2FILE": MOL2FILEPARM,
#     "CHARMMPSF": CHARMMPSF,
#     "CIFFILE": CIFFILE,
#     "GMXTOP": GMXTOP,
#     "SDFFILE": SDFFILE,
#     "TINKER": TINKERPARM,
#     "UNKNOWN_PARM": UNKNOWN_PARM,
# }

# Get the pytraj format key for the write_parm function for a specific file according to its format
def get_pytraj_parm_format (filename : str) -> str:
    if is_prmtop(filename):
        return 'AMBERPARM'
    if is_psf(filename):
        return 'CHARMMPSF'
    if is_top(filename):
        return 'GMXTOP'
    if is_pdb(filename):
        return 'PDBFILE'
    raise ValueError('The file ' + filename + ' format is not supported')