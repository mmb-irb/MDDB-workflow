
class NucleicAcid:
    """Attributes and methods for operations on nucleic acid sequences.

    :param sequence: nucleic acid sequence
    :type sequence: str
    :param unit_len: length of sequence subunit to be analysed
    :type unit_len: int
    :param unit_name: name of subunit to be analysed. If not given, it is inferred from `unit_len`.
    :type unit_name: str, optional
    :param ic_sequence: inverse-complement of nucleic acid sequence. If not given, it is computed from sequence.
    :type ic_sequence: str, optional
    :param Ac: complement to Adenine, defaults to "X"
    :type Ac: str, optional
    :raises ValueError: if nucleic acid sequence is empty
    :raises ValueError: if both U and T are present in sequence
    :raises ValueError: if unit length is negative or zero
    """

    def __init__(
            self,
            sequence,
            unit_len,
            unit_name="subunit",
            ic_sequence=None,
            Ac="X"):
        # sequence quality checks
        if unit_len <= 0:
            raise ValueError(
                "unit length must be an integer "
                "and greater than or equal to 1")
        if len(sequence) < unit_len:
            raise ValueError(
                "length of sequence can't be smaller than "
                "the specified subunit length")
        if "U" in sequence and "T" in sequence:
            raise ValueError("Both T and U bases present in sequence")
        self.sequence = sequence.upper()
        self.unit_len = unit_len
        self.unit_name = unit_name
        self._ic_sequence = ic_sequence
        self.Ac = self._adenine_complement(Ac)

        self.all_subunits = self.get_all_subunits()
        self.all_ic_subunits = self.get_all_ic_subunits()

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        representation = "NucleicAcid("
        for k, v in self.__dict__.items():
            representation += f"{k}={v}, "
        representation = representation[:-2] + ")"
        return representation

    @property
    def size(self):
        """get sequence length"""
        return len(self.sequence)

    @property
    def bases(self):
        """get nucleic acid sequence bases"""
        return set(self.sequence)

    @property
    def baselen(self):
        """get length of subunit base"""
        if self.unit_len % 2 == 0:
            baselen = 2
        else:
            baselen = 1
        return baselen

    @property
    def flanksize(self):
        """get size of subunit flanks"""
        return (self.unit_len - self.baselen) // 2

    @property
    def ic_sequence(self):
        """get inverse complement sequence"""
        if self._ic_sequence is None:
            self._ic_sequence = self.inverse_complement(self.sequence)
        return self._ic_sequence

    def _adenine_complement(self, Ac):
        """get complement base for Adenine (A)"""
        Ac = Ac.upper()
        if Ac != "T" and Ac != "U":
            # try to look for Adenine complement in sequence
            if "T" in self.sequence or "T" in self._ic_sequence:
                Ac = "T"
            elif "U" in self.sequence or "U" in self._ic_sequence:
                Ac = "U"
            else:
                Ac = "X"
                raise UserWarning(
                    "Adenine complement not U or T! Using X as placeholder "
                    "(this could lead to errors. "
                    "It is recommended to state this variable explicitly).")
        return Ac

    def inverse_complement(self, sequence):
        """compute inverse complement subunit."""
        complement = {
            "A": self.Ac,
            self.Ac: "A",
            "G": "C",
            "C": "G"
        }
        return "".join([complement[b] for b in sequence])[::-1]

    def get_subunit(self, idx):
        """
        Get complete subunit given the index of the first base, excluding flanks.
        For example, if you have the sequence AGTCTGA, and you want the tetramer at index 3, it will be TCTG.
        """
        subunit = self.sequence[
            idx - self.flanksize:
            idx + self.baselen + self.flanksize]
        return subunit

    def get_all_subunits(self):
        """
        Get a list of all subunits of length `unit_len` contained sequence.
        Returns a list of Subunit instances.
        """
        sequence_range = range(
            self.flanksize,
            self.size - self.baselen - self.flanksize + 1)
        return [self.get_subunit(idx) for idx in sequence_range]

    def get_all_ic_subunits(self):
        """Get a list of all subunits in sequence."""
        return [self.inverse_complement(unit) for unit in self.all_subunits]

    def search_subunit(self, subunit):
        """
        Returns index of the first base of subunit, excluding flanks (compatible with get_subunit method).
        For example, if you search for TCTG in sequence AGTCTGA, it will return index 3.
        """
        try:
            seq_index = self.sequence.index(subunit)
        except ValueError:
            try:
                seq_index = self.sequence.index(
                    self.inverse_complement(subunit))
            except ValueError:
                return None
        return seq_index + 1
