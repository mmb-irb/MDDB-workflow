from itertools import product, combinations_with_replacement
import numpy.ma as ma
import pandas as pd
import numpy as np
import os
from os import chdir, getcwd
from pathlib import Path

from mddb_workflow.tools.nassa_base import Base
from mddb_workflow.tools.nassa_loaders import load_sequence
from mddb_workflow.tools.nassa_loaders import load_serfile
from mddb_workflow.utils.constants import NASSA_ANALYSES_CANALS, BLUE_HEADER, COLOR_END, CYAN_HEADER, GREEN_HEADER

from mddb_workflow.utils.heatmaps_nassa import basepair_plot
from mddb_workflow.utils.heatmaps_nassa import bconf_heatmap
from mddb_workflow.utils.heatmaps_nassa import correlation_plot
from mddb_workflow.utils.heatmaps_nassa import arlequin_plot
from mddb_workflow.utils.bibitransformer_nassa import BiBiTransformer
from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.nassa_file import generate_nassa_config
from typing import Optional
import yaml


# There are five analyses in total

# BASE PAIR CORRELATION ANALYSIS
class BasePairCorrelation(Base):
    """Execution plan and methods for basepair correlation analysis."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def extract(self):
        extracted = {}
        sequences = []
        for seq_file in self.sequence_files:
            seq = load_sequence(
                seq_file,
                unit_name=self.unit_name,
                unit_len=self.unit_len)
            sequences.append(seq)
        self.Ac = seq.Ac
        self.logger.debug(f"Adenine complement set to: <{self.Ac}>")
        extracted["sequences"] = sequences
        self.logger.info(f"loaded {len(sequences)} sequences")

        #Iterate over the different parameters: shift, twist, roll, slide, tilt, rise
        bpcorr_coordiantes = NASSA_ANALYSES_CANALS['bpcorr']
        for bpcorr_coord in bpcorr_coordiantes:
            for coord in self.coordinate_info.keys():
                if coord == bpcorr_coord:
                    crd_data = []
                    for ser_file in self.coordinate_info[coord]:
                        ser = load_serfile(
                            ser_file,
                            self.tail,
                            self.n_lines)
                        crd_data.append(ser)
                    extracted[coord.lower()] = crd_data
                    self.logger.info(f"loaded {len(crd_data)} files for coordinate <{coord}>")

        return extracted

    def extraer_tetramero_central(hexamero):
                return hexamero[1:5]

    def transform(self, data):
        sequences = data.pop("sequences")
        # iterate over trajectories
        corr_results = {}
        for traj, seq in enumerate(sequences):
            trajectory_series = {coord.lower(): data[coord][traj]
                                 for coord in data.keys()}
            coordinate_corr = self.iterate_trajectory(
                seq, trajectory_series)
            corr_results[seq.sequence] = coordinate_corr
        joined_df = []
        for seq, val in corr_results.items():
            df = pd.DataFrame.from_dict(val).T
            joined_df.append(df)
        joined_df = pd.concat(joined_df)
        return joined_df

    def iterate_trajectory(self, sequence, coordinates):
        coordinate_correlations = {}
        # iterate over subunits
        start = 2 + sequence.flanksize
        end = sequence.size - (2 + sequence.baselen + sequence.flanksize - 1)
        # subtract 1 from end in order to stop at last sub1unit
        for idx in range(start, end-1):
            subunit = sequence.get_subunit(idx)
            next_subunit = sequence.get_subunit(idx+1)
            if subunit in coordinate_correlations:
                self.logger.info(
                    f"skipping repeated {self.unit_name} {subunit}...")
                continue
            self.logger.info(
                f"analyzing {self.unit_name} {subunit}...")
            # add 1 to idx since .ser table includes an index
            unit_df = pd.DataFrame(
                {coord: coordinates[coord][idx+1]
                 for coord in coordinates.keys()})
            next_unit_df = pd.DataFrame(
                {coord: coordinates[coord][idx+2]
                 for coord in coordinates.keys()})
            crd_corr = self.get_correlation(next_unit_df, unit_df)  # inverse order ?
            coordinate_correlations[f"{subunit}/{next_subunit}"] = crd_corr
        return coordinate_correlations

    def get_correlation(self, unit, next_unit):
        method = {
            "shift": "linear",
            "slide": "linear",
            "rise": "linear",
            "tilt": "circular",
            "roll": "circular",
            "twist": "circular"}
        coordinates = method.keys()
        combos = product(coordinates, repeat=2)
        result = {}
        for crd1, crd2 in combos:
            method1 = method[crd1]
            method2 = method[crd2]
            arr1 = unit[crd1]
            arr2 = next_unit[crd2]
            value = self.get_corr_by_method(
                method1,
                method2,
                arr1,
                arr2)
            result[f"{crd1}/{crd2}"] = value
        return result

    def get_corr_by_method(self, method1, method2, arr1, arr2):
        if method1 == "circular" and method2 == "linear":
            value = self.circlineal(arr2, arr1)
        if method1 == "linear" and method2 == "circular":
            value = self.circlineal(arr1, arr2)
        elif method1 == "circular" and method2 == "circular":
            value = self.circular(arr1, arr2)
        else:
            value = ma.corrcoef(ma.masked_invalid(arr1), ma.masked_invalid(arr2))[1, 0]
        return value

    @staticmethod
    def circular(x1, x2):
        x1 = x1 * np.pi / 180
        x2 = x2 * np.pi / 180
        diff_1 = np.sin(x1 - x1.mean())
        diff_2 = np.sin(x2 - x2.mean())
        num = (diff_1 * diff_2).sum()
        den = np.sqrt((diff_1 ** 2).sum() * (diff_2 ** 2).sum())
        return num / den

    @staticmethod
    def circlineal(x1, x2):
        x2 = x2 * np.pi / 180
        x1 = x1.to_numpy()
        x2 = x2.to_numpy()
        rc = ma.corrcoef(ma.masked_invalid(x1), ma.masked_invalid(np.cos(x2)))[1, 0]
        rs = ma.corrcoef(ma.masked_invalid(x1), ma.masked_invalid(np.sin(x2)))[1, 0]
        rcs = ma.corrcoef(ma.masked_invalid(np.sin(x2)), ma.masked_invalid(np.cos(x2)))[1, 0]
        num = (rc ** 2) + (rs ** 2) - 2 * rc * rs * rcs
        den = 1 - (rcs ** 2)
        correlation = np.sqrt(num / den)
        if ma.corrcoef(ma.masked_invalid(x1), ma.masked_invalid(x2))[1, 0] < 0:
            correlation *= -1
        return correlation

    def make_tables(self, dataset):
        dataset.to_csv(self.save_path / "all_basepairs.csv")

    def make_plots(self, dataset):
        basepair_plot(dataset, "all_basepairs", self.save_path)

# BASE PAIR CONFIRMATION ANALYSIS


class BConformations(Base):
    """Execution plan and methods for BI/BII conformations analysis pipeline"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        valid_coordinates = ["epsilC", "zetaC", "epsilW", "zetaW"]
        for coord in self.coordinate_info.keys():
            try:
                assert coord in valid_coordinates
            except AssertionError as e:
                raise ValueError(
                    f"{coord} is not a valid coordinate! "
                    "Please rename coordinates in your configuration file "
                    f"to match any of {valid_coordinates}") from e

    def extract(self):
        extracted = {}
        sequences = []
        for seq_file in self.sequence_files:
            seq = load_sequence(
                seq_file,
                unit_name=self.unit_name,
                unit_len=self.unit_len)
            sequences.append(seq)
        self.Ac = seq.Ac
        self.logger.debug(f"Adenine complement set to: <{self.Ac}>")
        extracted["sequences"] = sequences
        self.logger.info(f"loaded {len(sequences)} sequences")
        bconf_coordiantes = NASSA_ANALYSES_CANALS['bconf']
        for bconf_coord in bconf_coordiantes:
            for coord in self.coordinate_info.keys():
                if coord == bconf_coord:
                    crd_data = []
                    for ser_file in self.coordinate_info[coord]:
                        ser = load_serfile(
                            ser_file,
                            self.tail,
                            self.n_lines)
                        crd_data.append(ser)
                    extracted[coord.lower()] = crd_data
                    self.logger.info(
                        f"loaded {len(crd_data)} files for coordinate <{coord}>")
        return extracted

    def transform(self, data):
        sequences = data.pop("sequences")
        angles_df = []
        # get dataframe for each coordinate
        for traj, seq in enumerate(sequences):
            # start reading from the 4th column:
            # skip index, first two bases and first flanks
            start = 3 + seq.flanksize
            # skip last two bases
            end = seq.size - (seq.baselen + seq.flanksize)
            # select relevant subset of columns
            epsilC = data["epsilc"][traj][start:end]
            zetaC = data["zetac"][traj][start:end]
            epsilW = data["epsilw"][traj][start:end]
            zetaW = data["zetaw"][traj][start:end]
            traj_df = self.get_angles_difference(
                seq,
                epsilC,
                zetaC,
                epsilW,
                zetaW)
            angles_df.append(traj_df)
        angles_df = pd.concat(angles_df, axis=1)
        # percentages BI
        B_I = (angles_df < 0).sum(axis=0) * 100 / len(angles_df)  # self.n_lines
        # clean dataset
        B_I = B_I[~B_I.index.duplicated(keep='first')]
        B_I = B_I.reset_index()
        B_I = B_I.rename(columns={"index": self.unit_name, 0: "pct"})
        return B_I

    def get_angles_difference(self, seq, epsilC, zetaC, epsilW, zetaW):
        # get list of tetramers, except first and last two bases
        all_subunits = seq.all_subunits[2:-2]
        all_ic_subunits = seq.all_ic_subunits[2:-2]

        # concatenate zeta and epsil arrays
        zeta = pd.concat([
            pd.DataFrame(zetaW.T),
            pd.DataFrame(zetaC[::-1].T)],
            axis=1)
        zeta.columns = all_subunits * 2
        epsil = pd.concat([
            pd.DataFrame(epsilW.T),
            pd.DataFrame(epsilC[::-1].T)],
            axis=1)

        # This function was contained in a separated script in the original NASSA code (angle_utils.py)
        def fix_angle_range(x, domain=[0, 360]):
            """Fix angle range so it's in the given angle range (degrees) by adding or subtracting 360.

            :param float x: angle value (asumed to be in degrees)
            :param sequence domain: start and end of angle range.

            : return float: angle value with fixed range
            """
            while x < domain[0]:
                x += 360
            while x > domain[1]:
                x -= 360
            return x

        epsil.columns = all_subunits * 2
        # difference between epsilon and zeta coordinates
        diff = epsil - zeta
        diff = diff.applymap(lambda x: fix_angle_range(x, domain=[-180, 180]))

        # repeat with inverse-complementary sequence
        zeta_ic = pd.concat([
            pd.DataFrame(zetaW.T),
            pd.DataFrame(zetaC[::-1].T)],
            axis=1)
        zeta_ic.columns = all_ic_subunits * 2
        epsil_ic = pd.concat([
            pd.DataFrame(epsilW.T),
            pd.DataFrame(epsilC[::-1].T)],
            axis=1)
        epsil_ic.columns = all_ic_subunits * 2
        diff_ic = epsil_ic - zeta_ic
        diff_ic = diff_ic.applymap(
            lambda x: fix_angle_range(x, domain=[-180, 180]))

        diff_df = pd.concat([diff, diff_ic], axis=1)

        return diff_df

    def make_tables(self, B_I):
        B_I.to_csv(self.save_path / "BI.csv", index=False)
        # create BII
        B_II = B_I.copy()
        B_II["pct"] = 100 - B_I["pct"]
        B_II.to_csv(self.save_path / "BII.csv", index=False)

    def make_plots(self, B_I):
        bconf_heatmap(B_I, "BI", self.save_path, self.unit_len, self.Ac)
        B_II = B_I.copy()
        B_II["pct"] = 100 - B_I["pct"]
        bconf_heatmap(B_II, "BII", self.save_path, self.unit_len, self.Ac)


# COORDINATE CORRELATION ANALYSIS
class CoordinateCorrelation(Base):
    """Execution plan and methods for correlation analyses
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def extract(self):
        extracted = {}
        sequences = []
        for seq_file in self.sequence_files:
            seq = load_sequence(
                seq_file,
                unit_name=self.unit_name,
                unit_len=self.unit_len)
            sequences.append(seq)
        self.Ac = seq.Ac
        self.logger.debug(f"Adenine complement set to: <{self.Ac}>")
        extracted["sequences"] = sequences
        self.logger.info(f"loaded {len(sequences)} sequences")

        crdcorr_coordiantes = NASSA_ANALYSES_CANALS['crdcorr']
        for crdcorr_coord in crdcorr_coordiantes:
            for coord in self.coordinate_info.keys():
                if coord == crdcorr_coord:
                    crd_data = []
                    for ser_file in self.coordinate_info[coord]:
                        ser = load_serfile(
                            ser_file,
                            self.tail,
                            self.n_lines)
                        crd_data.append(ser)
                    extracted[coord.lower()] = crd_data
                    self.logger.info(
                        f"loaded {len(crd_data)} files for coordinate <{coord}>")
        return extracted

    def transform(self, data):
        sequences = data.pop("sequences")
        # iterate over trajectories
        corr_results = {}
        for traj, seq in enumerate(sequences):
            trajectory_series = {coord.lower(): data[coord][traj]
                                 for coord in data.keys()}
            correlations = self.iterate_trajectory(
                seq, trajectory_series)
            corr_results[traj] = correlations
        return corr_results

    def iterate_trajectory(self, sequence, coordinates):
        corrtype = {
            "shift": "linear",
            "slide": "linear",
            "rise": "linear",
            "tilt": "circular",
            "roll": "circular",
            "twist": "circular"}
        start = 3 + sequence.flanksize
        end = sequence.size - sequence.flanksize - sequence.baselen
        all_subunits = sequence.all_subunits[2:-2]
        results_df = {}
        for crd1, crd2 in combinations_with_replacement(
                coordinates.keys(),
                r=2):
            crd1_series = pd.DataFrame(data=coordinates[crd1][start:end]).T
            crd2_series = pd.DataFrame(data=coordinates[crd2][start:end]).T
            crd1_series.columns = all_subunits
            crd2_series.columns = all_subunits
            method = self.get_corr_method(corrtype[crd1], corrtype[crd2])
            corr_matrix = pd.DataFrame({
                col: crd1_series.corrwith(
                    crd2_series[col],
                    method=method
                ) for col in crd2_series.columns})
            crd1_inner_dict = results_df.get(crd1, {})
            crd1_inner_dict[crd2] = corr_matrix
            results_df[crd1] = crd1_inner_dict

        # build complete dataset
        coords = list(results_df.keys())
        for crd1, inner_dict in results_df.items():
            missing_coords = [
                crd for crd in coords if crd not in inner_dict.keys()]
            for crd2 in missing_coords:
                results_df[crd1][crd2] = results_df[crd2][crd1]
        dfs = {k: pd.concat(v) for k, v in results_df.items()}
        complete_df = pd.concat(dfs, axis=1)
        results_df["complete"] = complete_df
        return results_df

    def get_corr_method(self, corrtype1, corrtype2):
        if corrtype1 == "circular" and corrtype2 == "linear":
            method = self.circlineal
        if corrtype1 == "linear" and corrtype2 == "circular":
            method = self.circlineal
        elif corrtype1 == "circular" and corrtype2 == "circular":
            method = self.circular
        else:
            method = "pearson"
        return method

    @staticmethod
    def circular(x1, x2):
        x1 = x1 * np.pi / 180
        x2 = x2 * np.pi / 180
        diff_1 = np.sin(x1 - x1.mean())
        diff_2 = np.sin(x2 - x2.mean())
        num = (diff_1 * diff_2).sum()
        den = np.sqrt((diff_1 ** 2).sum() * (diff_2 ** 2).sum())
        return num / den

    @staticmethod
    def circlineal(x1, x2):
        x2 = x2 * np.pi / 180
        rc = np.corrcoef(x1, np.cos(x2))[1, 0]
        rs = np.corrcoef(x1, np.sin(x2))[1, 0]
        rcs = np.corrcoef(np.sin(x2), np.cos(x2))[1, 0]
        num = (rc ** 2) + (rs ** 2) - 2 * rc * rs * rcs
        den = 1 - (rcs ** 2)
        correlation = np.sqrt(num / den)
        if np.corrcoef(x1, x2)[1, 0] < 0:
            correlation *= -1
        return correlation

    def make_tables(self, dataset):
        for traj, data in dataset.items():
            # make a new directory to save data for each trajectory
            save_subpath = self.save_path / f"traj_{traj}"
            save_subpath.mkdir(exist_ok=True)
            # data.to_csv(save_subpath / "coordinate_corr.csv")
            # coords = data.columns.get_level_values(0).unique()
            for crd1, inner_dict in data.items():
                if crd1 == "complete":
                    inner_dict.to_csv(save_subpath / "coordinate_corr.csv")
                else:
                    for crd2, df in inner_dict.items():
                        df.to_csv(save_subpath / f"{crd1}_{crd2}.csv")

    def make_plots(self, dataset):
        complete_dataset = pd.concat(
            [value['complete'] for value in dataset.values()],
            axis=1)
        save_subpath = self.save_path / "complete"
        save_subpath.mkdir(exist_ok=True)
        correlation_plot(
            complete_dataset,
            "coordcorr",
            save_subpath)


# COORDINATE DISTRIBUTION ANALYSIS
class CoordinateDistributions(Base):
    """Execution plan and methods for coordinate distributions analysis pipeline."""

    def __init__(
            self,
            max_iter=400,
            tol=1e-5,
            **kwargs):
        super().__init__(**kwargs)
        self.max_iter = max_iter
        self.tol = tol
        self.logger.debug(f"max_iter: {self.max_iter}")
        self.logger.debug(f"tol: {self.tol}")

    def extract(self):
        extracted = {}
        sequences = []
        for seq_file in self.sequence_files:
            seq = load_sequence(
                seq_file,
                unit_name=self.unit_name,
                unit_len=self.unit_len)
            sequences.append(seq)
        self.Ac = seq.Ac
        self.logger.debug(f"Adenine complement set to: <{self.Ac}>")
        extracted["sequences"] = sequences
        self.logger.info(f"loaded {len(sequences)} sequences")

        coordist_coordiantes = NASSA_ANALYSES_CANALS['coordist']
        for coordist_coord in coordist_coordiantes:
            for coord in self.coordinate_info.keys():
                if coord == coordist_coord:
                    crd_data = []
                    for ser_file in self.coordinate_info[coord]:
                        ser = load_serfile(
                            ser_file,
                            self.tail,
                            self.n_lines)
                        crd_data.append(ser)
                    extracted[coord.lower()] = crd_data
                    self.logger.info(
                        f"loaded {len(crd_data)} files for coordinate <{coord}>")
        return extracted

    def transform(self, data):
        sequences = data.pop("sequences")
        subunits_classification = {}
        # get dataframe for each coordinate
        for coordinate, coord_dataset in data.items():
            coordinate_df = self.coordinate_iteration(
                sequences,
                coordinate,
                coord_dataset)
            coordinate_df = coordinate_df.drop_duplicates(
                subset=[self.unit_name])
            coordinate_df, global_mean, global_std = self.add_modality_labels(
                coordinate_df)
            subunits_classification[coordinate] = [
                coordinate_df, global_mean, global_std]
        return subunits_classification

    def coordinate_iteration(
            self,
            sequences,
            coordinate,
            coord_dataset):
        coordinate_df = []
        if self.duplicates:
            trajectory_df = self.duplicate_trajectory_iteration(
                coord_dataset,
                sequences,
                coordinate)
            coordinate_df.append(trajectory_df)
        else:
            for seq, dataseries in zip(sequences, coord_dataset):
                    trajectory_df = self.trajectory_iteration(
                        dataseries,
                        seq,
                        coordinate)
                    coordinate_df.append(trajectory_df)
        # concatenate all dataframes
        coordinate_df = pd.concat(
            coordinate_df,
            axis=0)
        return coordinate_df

    def duplicate_trajectory_iteration(
        self,
        coord_dataset,
        sequences,
        coordinate
    ):
        trajectory_info = []
        dict_all_subunits_ser = {}
        # iterate over datasets and sequences
        for dataseries, sequence in zip(coord_dataset, sequences):
            # iterate over subunits
            start = 2 + sequence.flanksize
            end = sequence.size - (2 + sequence.baselen + sequence.flanksize - 1)

            for idx in range(start, end):
                subunit = sequence.get_subunit(idx)
                ic_subunit = sequence.inverse_complement(subunit)
                ser = dataseries[idx + 1]
                ser = ser[~np.isnan(ser)].to_numpy()

                if subunit in dict_all_subunits_ser:
                    dict_all_subunits_ser[subunit] = np.concatenate((dict_all_subunits_ser[subunit], ser))
                else:
                    dict_all_subunits_ser[subunit] = ser

                if ic_subunit not in dict_all_subunits_ser:
                    dict_all_subunits_ser[ic_subunit] = ser

        for subunit, ser in dict_all_subunits_ser.items():
            ser = ser.reshape(-1, 1)
            subunit_information = self.subunit_iteration(
                ser,
                subunit,
                coordinate)
            if not self.bimod:
                subunit_information["unimodal"] = True
            trajectory_info.append(subunit_information)

        # create dataframe from list of dictionaries
        trajectory_df = pd.DataFrame.from_dict(trajectory_info)
        return trajectory_df

    def trajectory_iteration(
            self,
            dataseries,
            sequence,
            coordinate):
        trajectory_info = []
        # iterate over subunits
        start = 2 + sequence.flanksize
        end = sequence.size - (2 + sequence.baselen + sequence.flanksize - 1)
        for idx in range(start, end):
            # get unit and inverse-complement unit
            subunit = sequence.get_subunit(idx)
            ic_subunit = sequence.inverse_complement(subunit)
            self.logger.info(
                f"analyzing {self.unit_name} {subunit}/{ic_subunit}...")
            # add 1 to idx since .ser table includes an index
            ser = dataseries[idx + 1]
            ser = ser[~np.isnan(ser)].to_numpy()
            if ser.shape[0] < 2:
                self.logger.info("skipping because of insufficient data!")
                subunit_information = dict(
                    coordinate=coordinate,
                    binormal=False,
                    uninormal=True,
                    insuf_ev=False,
                    unimodal=True,
                    bics=[np.nan, np.nan],
                    mean1=np.nan,
                    mean2=np.nan,
                    var1=np.nan,
                    var2=np.nan,
                    w1=np.nan,
                    w2=np.nan)
                subunit_information[self.unit_name] = subunit
            else:
                # reshape dataset
                ser = ser.reshape(-1, 1)
                subunit_information = self.subunit_iteration(
                    ser,
                    subunit,
                    coordinate)
            if not self.bimod:
                subunit_information["unimodal"] = True
            trajectory_info.append(subunit_information)
            # add inverse-complement
            ic_subunit_information = subunit_information.copy()
            ic_subunit_information[self.unit_name] = ic_subunit
            trajectory_info.append(ic_subunit_information)
        # create dataframe from list of dictionaries
        trajectory_df = pd.DataFrame.from_dict(trajectory_info)
        return trajectory_df

    def subunit_iteration(
            self,
            ser,
            subunit,
            coordinate):
        # model subunit's series
        bibi = BiBiTransformer(max_iter=self.max_iter, tol=self.tol)
        subunit_info = bibi.fit_transform(ser)
        # add subunit name and coordinate to info dictionary
        subunit_info[self.unit_name] = subunit
        subunit_info["coordinate"] = coordinate
        return subunit_info

    def add_modality_labels(self, data):
        df = data.set_index(self.unit_name)

        # Series with global mean value for each tetramer
        weighted_mean = df.apply(
            lambda t: t["mean1"] if (
                t["uninormal"] or t["insuf_ev"] or np.isnan(t["mean2"])
            ) else (t["mean1"]*t["w1"]+t["mean2"]*t["w2"]),
            axis=1)

        # Series with comparison values
        comparison_1 = df.apply(
            lambda t: t["mean1"] if (
                t["uninormal"] or t["insuf_ev"] or not t["unimodal"]
            ) else t["mean1"]*t["w1"]+t["mean2"]*t["w2"],
            axis=1)

        # mean1: uninormal+unimodal
        # u1*w1+u2*w2: binormal+unimodal
        # mean2: binormal+bimodal
        comparison_2 = df.apply(
            lambda t: np.nan if t["unimodal"] else t["mean2"],
            axis=1)
        comparison_2 = comparison_2.fillna(comparison_1)

        # get limiting values for col1 and col2
        global_mean = weighted_mean.mean()
        global_std = weighted_mean.std()
        l1 = (global_mean + global_std)
        l2 = (global_mean - global_std)

        # 1: above mean+std (blue)
        # 0: between mean-std and mean+std (white)
        # -1: below mean-std (red)
        col1 = comparison_1.apply(
            lambda t: -1 if (t < l2) else (1 if (t > l1) else 0))
        col1 = col1.rename("col1")
        col2 = comparison_2.apply(
            lambda t: -1 if (t < l2) else (1 if (t > l1) else 0))
        col2 = col2.rename("col2")

        df = pd.concat(
            [df, col1, col2],
            axis=1).reset_index()

        return df, global_mean, global_std

    def make_tables(self, dataset, index=True):
        for coordinate, dataset in dataset.items():
            dataset[0].to_csv(
                self.save_path / f"{coordinate}.csv", index=index)

    def make_plots(self, dataset):
        for coordinate, data_dict in dataset.items():
            dataframes = data_dict[0]
            global_mean = data_dict[1]
            global_std = data_dict[2]
            arlequin_plot(
                dataframes,
                global_mean,
                global_std,
                coordinate,
                self.save_path,
                base=self.Ac,
                unit_name=self.unit_name,
                unit_len=self.unit_len)


# STIFFNESS ANALYSIS
class StiffnessDistributions(Base):

    def __init__(
            self,
            *args,
            **kwargs):
        super().__init__(
            *args,
            **kwargs)

    def extract(self):
        extracted = {}
        sequences = []
        for seq_file in self.sequence_files:
            seq = load_sequence(
                seq_file,
                unit_name=self.unit_name,
                unit_len=self.unit_len)
            sequences.append(seq)
        self.Ac = seq.Ac
        self.logger.debug(f"Adenine complement set to: <{self.Ac}>")
        extracted["sequences"] = sequences
        self.logger.info(f"loaded {len(sequences)} sequences")
        # In this case the selection of coordinates is different, it depends on the length of the unit (par or impar)
        if (self.unit_len % 2) == 0:
            stiff_coordiantes = NASSA_ANALYSES_CANALS['bpcorr']
        elif (self.unit_len % 2) == 1:
            stiff_coordiantes = NASSA_ANALYSES_CANALS['stiff']
        for stiff_coord in stiff_coordiantes:
            for coord in self.coordinate_info.keys():
                if coord == stiff_coord:
                    crd_data = []
                    for ser_file in self.coordinate_info[coord]:
                        ser = load_serfile(
                            ser_file,
                            self.tail,
                            self.n_lines)
                        crd_data.append(ser)
                    extracted[coord.lower()] = crd_data
                    self.logger.info(
                        f"loaded {len(crd_data)} files for coordinate <{coord}>")
        return extracted

    def transform(self, data):
        sequences = data.pop("sequences")
        results = {"stiffness": [], "covariances": {}, "constants": {}}
        for traj, seq in enumerate(sequences):
            traj_series = {coord.lower(): data[coord][traj]
                           for coord in data.keys()}
            traj_results = self.get_stiffness(
                seq,
                traj_series)
            results["stiffness"].append(traj_results["stiffness"])
            results["covariances"].update(traj_results["covariances"])
            results["constants"].update(traj_results["constants"])
        stiffness_df = pd.concat(results["stiffness"], axis=0)
        stiffness_df = stiffness_df.drop_duplicates(subset=[self.unit_name])
        stiffness_df = stiffness_df.set_index(self.unit_name)
        results["stiffness"] = stiffness_df
        return results

    def get_stiffness(
            self,
            sequence,
            series_dict):
        # get stiffness table for a given trajectory
        coordinates = list(series_dict.keys())
        results = {"stiffness": {}, "covariances": {}, "constants": {}}
        diagonals = {}
        start = 2 + sequence.flanksize
        end = sequence.size - (2 + sequence.baselen + sequence.flanksize - 1)
        for i in range(start, end):
            tetramer = sequence.get_subunit(i)
            ic_tetramer = sequence.inverse_complement(tetramer)
            cols_dict = {coord: series_dict[coord][i+1]
                         for coord in series_dict.keys()}
            stiffness_diag, cte, cov_df = self.get_subunit_stiffness(
                cols_dict,
                coordinates)
            diagonals[tetramer] = np.append(
                stiffness_diag,
                [np.prod(stiffness_diag), np.sum(stiffness_diag)])
            diagonals[ic_tetramer] = np.append(
                stiffness_diag,
                [np.prod(stiffness_diag), np.sum(stiffness_diag)])
            # results["covariances"][tetramer] = cov_df
            results["covariances"][ic_tetramer] = cov_df
            # results["constants"][tetramer] = cte
            results["constants"][ic_tetramer] = cte
        # build stiffness table
        columns = [sequence.unit_name] + coordinates + ["product", "sum"]
        results["stiffness"] = pd.DataFrame.from_dict(
            diagonals,
            orient="index").reset_index()
        results["stiffness"].columns = columns
        return results

    def get_subunit_stiffness(
            self,
            cols_dict,
            coordinates,
            scaling=[1, 1, 1, 10.6, 10.6, 10.6],
            KT=0.592186827):
        if (self.unit_len % 2) == 0:
            SH_av = cols_dict["shift"].mean()
            SL_av = cols_dict["slide"].mean()
            RS_av = cols_dict["rise"].mean()
            TL_av = self.circ_avg(cols_dict["tilt"])
            RL_av = self.circ_avg(cols_dict["roll"])
            TW_av = self.circ_avg(cols_dict["twist"])
        elif (self.unit_len % 2) == 1:
            SH_av = cols_dict["shear"].mean()
            SL_av = cols_dict["stretch"].mean()
            RS_av = cols_dict["stagger"].mean()
            CW_av = cols_dict["chiw"].mean()
            CC_av = cols_dict["chic"].mean()
            TL_av = self.circ_avg(cols_dict["buckle"])
            RL_av = self.circ_avg(cols_dict["propel"])
            TW_av = self.circ_avg(cols_dict["opening"])
        cols_arr = [cols_dict[coord] for coord in coordinates]
        cols_arr = np.array(cols_arr).T

        cv = ma.cov(ma.masked_invalid(cols_arr), rowvar=False)
        cv.filled(np.nan)

        cov_df = pd.DataFrame(cv, columns=coordinates, index=coordinates)
        stiff = np.linalg.inv(cv) * KT
        # Added two new variables: ChiC and ChiW -> 8 (for PENTAMERS)
        if (self.unit_len % 2) == 0:
            last_row = [SH_av, SL_av, RS_av, TL_av, RL_av, TW_av]  #, CW_av, CC_av]
            stiff = np.append(stiff, last_row).reshape(7, 6)
        elif (self.unit_len % 2) == 1:
            last_row = [SH_av, SL_av, RS_av, TL_av, RL_av, TW_av, CW_av, CC_av]
            stiff = np.append(stiff, last_row).reshape(9, 8)
            scaling = [1, 1, 1, 10.6, 10.6, 10.6, 1, 1]

        stiff = stiff.round(6)
        stiff_diag = np.diagonal(stiff) * np.array(scaling)

        cte = pd.DataFrame(stiff)
        cte.columns = coordinates
        cte.index = coordinates + ["avg"]
        return stiff_diag, cte, cov_df



    @ staticmethod
    def circ_avg(xarr, degrees=True):
        n = len(xarr)
        if degrees:
            # convert to radians
            xarr = xarr * np.pi / 180
        av = np.arctan2(
            (np.sum(np.sin(xarr)))/n,
            (np.sum(np.cos(xarr)))/n) * 180 / np.pi
        return av

    def unimod_labels(self, df):
        # get limiting values for col1 and col2
        global_mean = df.mean(axis=0)
        global_std = df.std(axis=0)
        l1 = global_mean + global_std
        l2 = global_mean - global_std

        # 1: above mean+std (blue)
        # 0: between mean-std and mean+std (white)
        # -1: below mean-std (red)
        labeled_df = (df < l2) * -1 + (df > l1)
        labeled_df.loc["g_mean"] = global_mean
        labeled_df.loc["g_std"] = global_std

        return labeled_df

    def make_tables(self, dataset):
        # stiffness
        stiffness_data = dataset["stiffness"]
        stiffness_data.to_csv(self.save_path / "stiffness.csv")
        # covariances
        covariances_path = self.save_path / "covariances"
        covariances_path.mkdir(exist_ok=True)
        for key, val in dataset["covariances"].items():
            val.to_csv(covariances_path / f"{key}.csv")
        # constants
        constants_path = self.save_path / "constants"
        constants_path.mkdir(exist_ok=True)
        for key, val in dataset["constants"].items():
            val.to_csv(constants_path / f"{key}.csv")

    def make_plots(self, dataset):
        stiffness_data = dataset["stiffness"]
        labeled_df = self.unimod_labels(stiffness_data)
        for col in labeled_df.columns:
            df = labeled_df[col]
            g_mean = df.loc["g_mean"]
            g_std = df.loc["g_std"]
            df = df.iloc[:-2]
            df = df.rename("col1")
            df = df.reset_index()
            df["col2"] = df["col1"].copy()
            arlequin_plot(
                df,
                g_mean,
                g_std,
                col,
                self.save_path,
                base=self.Ac,
                unit_name=self.unit_name,
                unit_len=self.unit_len)


# Run the NASSA analysis pipeline giving the name of the analysis and the configuration file.
# The overwrite flag could be used to overwrite the output folder if it already exists.
def run_nassa(analysis_name: str,
              config_archive: dict,
              overwrite_nassa: bool = False):
    """Run the NASSA analysis pipeline."""
    # Dictionary with the available NASSA analyses and their corresponding classes
    analyses = {
        "bpcorr": BasePairCorrelation,
        "bconf": BConformations,
        "crdcorr": CoordinateCorrelation,
        "coordist": CoordinateDistributions,
        "stiff": StiffnessDistributions
    }

    # We must check if the analysis name is valid
    if analysis_name in analyses:
        # A folder for each analysis has to be created
        analyses_folder = os.path.join(config_archive["save_path"], analysis_name)
        if not os.path.exists(analyses_folder):
            os.mkdir(analyses_folder)
            # The output folder is added to the configuration archive
            config_archive["save_path"] = analyses_folder
        else:
            # If the output folder already exists, it is checked so that it can be overwritten with te overwrite_nassa flag
            if os.listdir(analyses_folder) == [] or overwrite_nassa:
                config_archive["save_path"] = analyses_folder
                print(f'WARNING: Output folder {analyses_folder} already exists. Overwriting it.')
            else:
                print(f'WARNING: Output folder {analyses_folder} already exists and is not empty. Skipping NASSA analysis. \nSet the overwrite flag to overwrite the output folder (-own).')
                return
        # The analysis is run
        analysis_class = analyses[analysis_name]
        # The stiffness analysis is a special case, since the coordinate files are different depending on the unit length (if it is pair or impair)
        if analysis_name == "stiff":
            unit_len = config_archive["unit_len"]  # Obtain the unit length from the configuration archive
            if (unit_len % 2) == 0:
                coordinate_files = NASSA_ANALYSES_CANALS["bpcorr"]
            elif (unit_len % 2) == 1:
                coordinate_files = NASSA_ANALYSES_CANALS["stiff"]
        # The rest of the analyses have their corresponding coordinate files
        else:
            coordinate_files = NASSA_ANALYSES_CANALS[analysis_name]
        # The coordinate files are filtered from the configuration archive because we only need the ones that correspond to the analysis
        config_archive["coordinate_info"] = {
            coord: config_archive["coordinate_info"][coord] for coord in coordinate_files}
        # Call the analysis class with the configuration archive and run NASSA software
        analysis_instance = analysis_class(**config_archive)
        analysis_instance.run()
    # If the analysis name is not valid, an error is raised
    else:
        raise ValueError(
            f"{analysis_name} is not a valid analysis! "
            "Please choose from: "
            f"{list(analyses.keys())}")


def workflow_nassa(
    config_file_path: Optional[str],
    analysis_names: Optional[list[str]],
    make_config: bool = False,
    output: Optional[list[str]] = None,
    working_directory: str = '.',
    overwrite: bool = False,
    overwrite_nassa: bool = False,
    helical_par: bool = False,
    accession: Optional[str] = None,
    sample_trajectory: bool = False,
    proj_dirs: Optional[list[str]] = None,
    input_structure_file: Optional[str] = None,
    input_trajectory_file: Optional[str] = None,
    input_top_file: Optional[str] = None,
    all: bool = False,
    unit_len: int = 6,
    n_sequences: Optional[int] = '*',
    seq_path: str = None,
    md_directories: Optional[list[str]] = None,
    trust: bool = False,
    mercy: bool = False,
    ):

    # Change to the working directory and print the name of the directory to inform the user
    chdir(working_directory)
    current_directory_name = getcwd().split('/')[-1]
    print(f'\n{BLUE_HEADER}Running NASSA at {current_directory_name}{COLOR_END}')
    # If the user select the helical parameter analysis, we need to call the function from workflow to run it
    if helical_par:
        # AGUS: This option with glob is commented because it is not working properly so if in the future it is needed, it should be fixed
        # AGUS: the objective to use glob is to run this part (of project directories) calling this function from python and not from the command line
        # seq_paths = glob.glob(f'{proj_dirs}')
        # seq_paths = [path for path in glob.glob(f'{proj_dirs}*') if os.path.isdir(path)]

        actual_path = getcwd()
        # The workflow function is imported from mddb_workflow to run the helical parameters analysis and do all the checks in each project directory
        # By this way, we are assuming that each project directory is an independent project or MD simulation
        from mddb_workflow.mwf import workflow
        # Iterate over the project directories
        for proj_path in proj_dirs:
            # Obtain the path of the project directory and change to it
            md_path = os.path.join(actual_path, proj_path)
            os.chdir(md_path)
            print(f'\n{CYAN_HEADER}Running Helical Parameters at {proj_path}{COLOR_END}')
            # Call the workflow function to run the helical parameters analysis with the include flag set to helical and overwrite (if it is added)
            workflow(project_parameters={
                'accession': accession,
                'sample_trajectory': sample_trajectory,
                'input_topology_filepath': input_top_file,
                'input_trajectory_filepaths': input_trajectory_file,
                'md_directories': md_directories,
                'trust': trust,
                'mercy': mercy
                },
                include=['helical'],
                overwrite=overwrite
            )
            os.chdir('..')
        # If all flag is selected, the NASSA analysis will be run after the helical parameters analysis
        # Reminder: the NASSA analysis needs the helical parameters files to be generated before running it.
        if all:
            # The sequences path is needed to run the NASSA analysis so it has to be defined
            if seq_path is None:
                raise InputError('No sequence path defined. Please define it as an argument (--seq_path)')
            output_nassa_analysis = os.path.join(actual_path, 'nassa_analysis')
            # If the output folder already exists, it is checked so that it can be overwritten with te overwrite_nassa flag
            if not os.path.exists(output_nassa_analysis):
                os.mkdir(output_nassa_analysis)
            else:
                if overwrite_nassa is True:
                    print(f'WARNING: Output folder {output_nassa_analysis} already exists. Overwriting it.')
                else:
                    print(f'WARNING: Output folder {output_nassa_analysis} already exists and is not empty. Skipping NASSA analysis. \nSet the overwrite flag to overwrite the output folder (--overwrite_nassa).')
                    return
            # Generate the NASSA configuration file
            generate_nassa_config(
                    folder_path=proj_dirs,
                    seq_path=seq_path,
                    output_path=output_nassa_analysis,
                    unit_len=unit_len,
                    n_sequences=n_sequences,
            )
            print(f'NASSA configuration file generated at {output_nassa_analysis}')
            config_file_path = os.path.join(output_nassa_analysis, 'nassa.yml')
            # The user can select the analysis to run
            if analysis_names is not None:
                for analysis_name in analysis_names:
                    print(f'  Running analysis {analysis_name}')
                    # Read the configuration file
                    try:
                        with Path(config_file_path).open("r") as ymlfile:
                            config_archive = yaml.load(
                                ymlfile, Loader=yaml.FullLoader)
                    except FileNotFoundError:
                        raise FileNotFoundError(
                            f"Configuration file {config_file_path} not found!")
                    # Run the NASSA analysis with the selected analysis
                    run_nassa(analysis_name=analysis_name,
                            config_archive=config_archive,
                            overwrite_nassa=overwrite_nassa)
                print(f'NASSA analysis completed at {output_nassa_analysis}')
                return
            # If the user does not select any analysis, all the analyses will be run
            for analysis_name in NASSA_ANALYSES_CANALS.keys():
                print(f'  Running analysis {analysis_name}')
                # Read the configuration file
                try:
                    with Path(config_file_path).open("r") as ymlfile:
                        config_archive = yaml.load(
                            ymlfile, Loader=yaml.FullLoader)
                except FileNotFoundError:
                    raise FileNotFoundError(
                        f"Configuration file {config_file_path} not found!")
                # Run the NASSA analysis with the selected analysis
                run_nassa(analysis_name=analysis_name,
                          config_archive=config_archive,
                          overwrite_nassa=overwrite_nassa)
            print(f'NASSA analysis completed at {output_nassa_analysis}')
            return
    # If the helical parameter analysis is not selected, the NASSA analysis will be run
    # To run this, the user has to provide the configuration file with the information needed
    if config_file_path and analysis_names is not None:
        config_file_path = Path(config_file_path)
        print(f'  Using config file {config_file_path}')
        # Load the configuration file
        try:
            with Path(config_file_path).open("r") as ymlfile:
                config_archive = yaml.load(
                    ymlfile, Loader=yaml.FullLoader)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Configuration file {config_file_path} not found!")
        # The user can select the analysis to run
        for analysis_name in analysis_names:
            # THe output folder is defined in the configuration file
            if output is not None:
                output_path = os.path.join(output, analysis_name)
                config_archive['save_path'] = output_path
            else:
                if 'save_path' in config_archive:
                    output_path = os.path.join(config_archive['save_path'], analysis_name)
                    config_archive['save_path'] = output_path
                else:
                    raise InputError('No output path defined. Please define it in the configuration file or pass it as an argument (--output)')

            # To check if the user has created the configuration file with the needed information
            # it's important to check if the needed files are defined for the analysis selected
            canal_files_config = config_archive['coordinate_info']
            needed_files = NASSA_ANALYSES_CANALS[analysis_name]
            for file in needed_files:
                # If some of the files needed are not in the config file, raise an error because it is needed
                if file not in canal_files_config:
                    raise InputError(f'Analysis {analysis_name} requires the files of {file} coordinate to be defined in the configuration file')

            # The check for the output folder is done here
            if os.path.exists(output_path):
                if os.listdir(output_path) == [] or overwrite_nassa == True:
                    print(f'  Output folder {output_path} already exists. Overwriting analysis {analysis_name}')
                    run_nassa(analysis_name, config_archive, overwrite_nassa)
                else:
                    print(f'  Output folder {output_path} already exists. Skipping analysis. \nSet the overwrite flag to overwrite the output folder (--overwrite_nassa).')
                    continue
            else:
                os.mkdir(output_path)
                print(f'  Running analysis {analysis_name} and saving results in {output_path}')
                run_nassa(analysis_name, config_archive, overwrite_nassa)
    if all:
        for analysis_name in NASSA_ANALYSES_CANALS.keys():
            print(f"{GREEN_HEADER} |-----> Running analysis {analysis_name}{COLOR_END}")
            # Read the configuration file
            try:
                with Path(config_file_path).open("r") as ymlfile:
                    config_archive = yaml.load(
                        ymlfile, Loader=yaml.FullLoader)
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Configuration file {config_file_path} not found!")
            # Override nassa.yml output path if defined in the command line
            if output is not None:
                output_path = os.path.join(output, 'nassa_analysis')
            else:
                output_path = os.path.join(config_archive['save_path'], 'nassa_analysis')
            config_archive['save_path'] = output_path
            # Check if the NASSA output folder exists
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            # Check if the analysis output folder exists in NASSA output folder
            if not os.path.exists(os.path.join(output_path, analysis_name)):
                os.makedirs(os.path.join(output_path, analysis_name))
            #     config_archive['save_path'] = os.path.join(output_path, analysis_name)
            # Run the NASSA analysis with the selected analysis
            run_nassa(analysis_name=analysis_name,
                        config_archive=config_archive,
                        overwrite_nassa=overwrite_nassa)
            print(f'NASSA analysis completed at {current_directory_name}')
