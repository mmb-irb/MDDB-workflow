import pathlib
from collections import deque
from mddb_workflow.utils.nucleicacid import NucleicAcid
import pandas as pd

def load_sequence(seqfile, unit_len, unit_name=None):
    """
    Load single text file containing forward and inverse-complementary sequences.

    :param str seq_file: a string with the path to the sequences file.
    :param str unit_name: name of subunit to be analyzed.
    :param int unit_len: length of subunit to be analyzed.
    :returns: NASequence object for loaded sequence.
    """
    assert isinstance(seqfile, str) or isinstance(seqfile, pathlib.Path)
    sequences = pathlib.Path(seqfile).read_text().split()
    if len(sequences) < 2:
        sequences.append(None)
        try:
            assert len(sequences) == 2
        except AssertionError:
            raise AssertionError("Error in sequence file! Check its not empty")
    nucleic_acid = NucleicAcid(
        sequence=sequences[0],
        ic_sequence=sequences[1],
        unit_name=unit_name,
        unit_len=unit_len)
    return nucleic_acid

# This function is a copy of the previous one, but it does not read from a file.
# It is used to load a sequence from a string. It is used in the helical parameters analysis.
def load_sequence2(sequence, unit_len, unit_name='hexamer'):
    """
    Load single text file containing forward and inverse-complementary sequences.

    :param str seq_file: a string with the path to the sequences file.
    :param str unit_name: name of subunit to be analyzed.
    :param int unit_len: length of subunit to be analyzed.
    :returns: NASequence object for loaded sequence.
    """
    #assert isinstance(seqfile, str) or isinstance(seqfile, pathlib.Path)
    #sequences = pathlib.Path(seqfile).read_text().split()
    sequences = [sequence]
    if len(sequences) == 1:
        sequences.append(reverse_sequence(sequence))
    if len(sequences) < 2:
        sequences.append(None)
        try:
            assert len(sequences) == 2
        except AssertionError:
            raise AssertionError("Error in sequence file! Check its not empty")
    nucleic_acid = NucleicAcid(
        sequence=sequences[0],
        ic_sequence=sequences[1],
        unit_name=unit_name,
        unit_len=unit_len)
    return nucleic_acid
# This function is used to reverse the sequence and obtain the inverse complement.
def reverse_sequence(sequence,DNA=True):
    if DNA: # If DNA flag is tru we want to compute the inverse using T instead of U
        A_base = "T"
    else: # Now it is RNA so we want to convert T to U
        A_base = "U"
    inverse = {"A":A_base,"G":"C","C":"G",A_base:"A"} # Dictionary to convert easily the sequence 
    inv_seq = ""
    for i in sequence[::-1]: # Traverse the sequence from the end to the beginning
        inv_seq += inverse[i] # Obtain the complementary base 
    return inv_seq # Return the inverse complement

def write_sequence(nucleic_acid, filename):
    assert isinstance(nucleic_acid, NucleicAcid)
    output = f"{nucleic_acid.sequence}\n{nucleic_acid.ic_sequence}"
    pathlib.Path(filename).write_text(output)



def load_serfile(ser_file, tail=True, n_lines=None):
    """
    Load single file containing a coordinate's series.

    :param str ser_file: path to .ser file.
    :param bool tail: (Default True) read the last ``n_lines`` of the file. Otherwise, read the first ``n_lines``.
    :param int n_lines: number of rows to read.
    :returns pandas.DataFrame: .ser file converted into a pandas.DataFrame table
    """
    if tail:
        #with open(ser_file, "r") as f:
            # read last line of file, get index number for that line
        #    total_lines = int(deque(f, 1)[0].split()[0])
        #extra_kwargs = dict(skiprows=total_lines - n_lines)
        with open(ser_file, "r") as f:
            total_lines_str = deque(f, 1)[0].split()[0]
            total_lines = int(total_lines_str) if total_lines_str.isdigit() else None

        if total_lines is not None and n_lines is not None:
            extra_kwargs = dict(skiprows=max(0, total_lines - n_lines))
        else:
            extra_kwargs = dict()
    else:
        extra_kwargs = dict(nrows=n_lines)
    ser_data = pd.read_csv(
        ser_file,
        header=None,
        sep='\s+',
        index_col=0,
        **extra_kwargs)
    return ser_data

def write_serfile(data, filename, indent=8, decimals=2, transpose=True):
    """Write data to same format as .ser file.
    By default, data is asumed to be in shape (n_cols, n_frames), and it's written in 8-spaced columns with values rounded to two decimals.

    :param numpy.ndarray data: output data
    :param str filename: dataset's filename
    :param indent: width of columns, defaults to 8
    :type indent: int, optional
    :param decimals: number of rounding decimals, defaults to 2
    :type decimals: int, optional
    :param transpose: transpose data array before writing. It should be used so array shape is (n_frames, n_cols). Defaults to True
    :type transpose: bool, optional
    """
    if transpose:
        data = data.T
    with open(filename, "w") as f:
        for row in data:
            s = f"{int(row[0]):>indent}"
            for elem in row[1:]:
                elem = round(elem, decimals)
                s += f"{elem:>indent}"
            s += "\n"
            f.write(s)
