# Protein residues
aminoacids = {
    "ALA":"A",
    "ALAN":"A",
    "ALAC":"A",
    "ARG":"R",
    "ARGN":"R",
    "ARGC":"R",
    "ASN":"N",
    "ASNN":"N",
    "ASNC":"N",
    "ASP":"D",
    "ASPN":"D",
    "ASPC":"D",
    "CYS":"C",
    "CYSN":"C",
    "CYSC":"C",
    "CYH":"C",
    "CSH":"C",
    "CSS":"C",
    "CYX":"C",
    "CYP":"C",
    "GLN":"Q",
    "GLNN":"Q",
    "GLNC":"Q",
    "GLU":"E",
    "GLUN":"E",
    "GLUC":"E",
    "GLY":"G",
    "GLYN":"G",
    "GLYC":"G",
    "HIS":"H",
    "HISN":"H",
    "HISC":"H",
    "HID":"H",
    "HIE":"H",
    "HIP":"H",
    "HSD":"H",
    "HSE":"H",
    "ILE":"I",
    "ILEN":"I",
    "ILEC":"I",
    "ILU":"I",
    "LEU":"L",
    "LEUN":"L",
    "LEUC":"L",
    "LYS":"K",
    "LYSN":"K",
    "LYSC":"K",
    "MET":"M",
    "METN":"M",
    "METC":"M",
    "PHE":"F",
    "PHEN":"F",
    "PHEC":"F",
    "PRO":"P",
    "PRON":"P",
    "PROC":"P",
    "PRØ":"P",
    "PR0":"P",
    "PRZ":"P",
    "SER":"S",
    "SERN":"S",
    "SERC":"S",
    "THR":"T",
    "THRN":"T",
    "THRC":"R",
    "TRP":"W",
    "TRPN":"W",
    "TRPC":"W",
    "TRY":"W",
    "TYR":"Y",
    "TYRN":"Y",
    "TYRC":"Y",
    "VAL":"V",
    "VALN":"V",
    "VALC":"V"
}

# Nucleic acid residues
nucleotides = {
    "A": "A",
    "A3": "A",
    "A5": "A",
    "C": "C",
    "C3": "C",
    "C5": "C",
    "T": "T",
    "T3": "T",
    "T5": "T",
    "G": "G",
    "G3": "G",
    "G5": "G",
    "U": "U",
    "U3": "U",
    "U5": "U",
    "DA": "A",
    "DT": "T",
    "DC": "C",
    "DG": "G",
}

# Transform a residue name to its equivalent single letter code
# If the residue name is not recognized then return "X"
# e.g. "ARG" -> "R", "WTF" -> "X"
# You can choose which residue types are targeted (e.g. aminoacids only)
# Options are: 'all', 'aminoacids' or 'nucleotides'
# All residue types are targeted by default
def residue_name_2_letter (residue_name : str, residue_types : str = "all") -> str:
    # Set the target residues
    if residue_types == "all":
        target_residues = { **aminoacids, **nucleotides }
    elif residue_types == "aminoacids":
        target_residues = aminoacids
    elif residue_types == "nucleotides":
        target_residues = nucleotides
    else:
        raise ValueError('Unrecognized residue types ' + str(residue_types))
    # Now find the corresponding letter among the target residues
    ref = target_residues.get(residue_name, False)
    return ref if ref else "X"