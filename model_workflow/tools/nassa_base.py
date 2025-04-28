import time
import pathlib
import logging
from abc import ABC, abstractmethod


class Base(ABC):
    """Base class for nucleic acid analysis workflow.

    :param list sequence_files: paths to sequence files.
    :param dict coordinate_info: dictionary with coordinates as keys, and coordinate series files as values.
    :param int n_lines: number of lines to read from coordinate series files.
    :param bool tail: read last ``n_lines`` from coordinate series files.
    :param str unit_name: name of subunit to analyze, used to name data table columns.
    :param int unit_len: length of subunit to analyze.
    :param bool save_tables: save data tables as .csv files.
    :param bool save_plots: save data visualizations as .pdf files.
    :param bool verbose: print verbose output.
    :param str save_path: path to save directory. If it doesn't exist, it is created during execution. Default is the current working directory.

    :raises ValueError: If number of sequence files doesn't match number of coordinate files.
    """

    def __init__(
            self,
            sequence_files,
            coordinate_info,
            save_path,
            unit_len,
            unit_name="subunit",
            n_lines=1000,
            tail=True,
            bimod=True,
            save_tables=True,
            save_plots=True,
            verbose=True):

        # sequence and coordinate paths and files from config file
        for coordinate_files in coordinate_info.values():
            if (len(sequence_files) != len(coordinate_files)):
                raise ValueError(
                    "number of sequence files must match "
                    "number of coordinate files")
        self.sequence_files = sequence_files
        self.coordinate_info = coordinate_info

        # variables
        self.n_lines = n_lines
        self.tail = tail
        self.unit_name = unit_name
        self.unit_len = unit_len
        self.bimod = bimod

        # create logger
        self.logger = self.create_logger(verbose)

        # parse paths from string to PosixPath
        self._save_path = save_path

        # flags
        self.save_tables = save_tables
        self.save_plots = save_plots

        # log information
        self.logger.info(
            "number of sequence files: "
            f"{len(self.sequence_files)}")
        self.logger.info(
            "number of files for each coordinate: "
            f"{[(k, len(v)) for k, v in self.coordinate_info.items()]}")
        self.logger.debug(f"sequence files: {self.sequence_files}")
        self.logger.debug(f"coordinate files: {self.coordinate_info}")
        self.logger.debug(
            f"reading {n_lines} lines from "
            f"{'tail' if self.tail else 'head'} of input files...")
        self.logger.debug(
            f"analyzing units of length {self.unit_len}...")

    @ abstractmethod
    def extract():
        """Extract data to be analyzed"""
        pass

    @ abstractmethod
    def transform():
        """Perform data transformations, cleaning and analyses"""
        pass

    @ abstractmethod
    def make_tables(data, **kwargs):
        """Save data in tables"""
        pass

    @ abstractmethod
    def make_plots(data, **kwargs):
        """Save data visualizations"""
        pass

    def load(self, data, **kwargs):
        """Save data in table and visualization formats.

        :param data: processed datasets
        """
        self.make_tables(data, **kwargs)
        self.make_plots(data, **kwargs)

    def run(self, **kwargs):
        """Run complete data analysis"""
        start = time.time()

        data = self.extract(**kwargs)
        data = self.transform(data, **kwargs)

        if self.save_tables:
            self.make_tables(data, **kwargs)
        if self.save_plots:
            self.make_plots(data, **kwargs)
        end = time.time()
        self.logger.info(
            f"Execution took {end-start:.2f} seconds "
            f"({(end-start)//60:.0f} minutes).")

        return data

    @ property
    def save_path(self):
        """
        Parse ``save_path`` directory. If it doesn't exist, directory is created along with parent directories.
        If not provided, current directory is used.

        :return pathlib.Path: Save directory path.
        """
        if self._save_path:
            new_path = pathlib.Path(self._save_path)
            if not new_path.exists():
                self.logger.info("creating directory to save output files...")
                self.logger.debug(f"created path {new_path}")
                new_path.mkdir(parents=True, exist_ok=True)
            return new_path
        else:
            if self.save_plots or self.save_tables:
                self.logger.info(
                    "save path not provided, using current directory")
                return pathlib.Path.cwd()
            else:
                return None

    @ staticmethod
    def create_logger(verbose):
        """Create logger.

        :param bool verbose: if True, logging level is DEBUG. Else, it's set to INFO.
        :return RootLogger: logger
        """
        if verbose:
            logging_level = logging.DEBUG
        else:
            logging_level = logging.INFO
        logging_format = "[%(filename)s:%(lineno)s] %(levelname)s: %(message)s"
        logging.basicConfig(level=logging_level, format=logging_format)
        logger = logging.getLogger('NASSA')
        return logger
