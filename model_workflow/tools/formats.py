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
    raise ValueError('The file ' + filename + ' format is not supported')