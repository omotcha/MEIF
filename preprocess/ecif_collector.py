"""
platform: any
env: any
name: ecif_collector.py
collect ECIF fingerprints of all protein-ligand pairs
"""
import os.path

import pandas as pd
from configs.config import data_dir, meif_data_dir, dataset_2016_dir, aff_data_file, data_aug_dir
from util.ECIF import ECIF, LIGAND_DESC
import csv


class ECIF_Collector:
    _distance_cutoffs = []
    _ecif_helper = None

    def __init__(self, dc_list):
        if dc_list is None or len(dc_list) == 0:
            self._distance_cutoffs = [6.0]
        else:
            self._distance_cutoffs = dc_list
        self._ecif_helper = ECIF()

    # callables

    def collect_ecif(self, tag):
        """
        collect ECIF data from pdb files and sdf files
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        :return:
        """
        def calc_ecif(iname, itag, distance_cutoff):
            if itag == "core":
                return
            if itag == "general":
                itag = "general-minus-refined"
            protein_file = os.path.join(dataset_2016_dir[itag], iname, "{}_protein.pdb".format(iname))
            ligand_file = os.path.join(dataset_2016_dir[itag], iname, "{}_ligand.sdf".format(iname))
            ecif_list = self._ecif_helper.get_ecif(protein_file,
                                                   ligand_file,
                                                   float(distance_cutoff))
            ecif_writer.writerow([iname] + ecif_list)
            return

        print("\nCollecting ECIF data: \n")
        aff_data = pd.read_csv(os.path.join(data_dir, "affinity_data_2016.csv"))
        # for functional testing
        # aff_data = aff_data[0:50]
        for d in self._distance_cutoffs:
            print("\n distance_cutoff: {}\n".format(d))
            if tag is None or tag == "":
                f_ecif = open(os.path.join(meif_data_dir, "ECIF_{}.csv".format(d)), 'a', newline='')
            else:
                f_ecif = open(os.path.join(meif_data_dir, "ECIF_{}_{}.csv".format(d, tag)), 'a', newline='')
            # write csv header
            ecif_header = ["PDB"] + self._ecif_helper.get_possible_pl()
            ecif_writer = csv.writer(f_ecif)
            ecif_writer.writerow(ecif_header)
            list(aff_data.apply(lambda x: calc_ecif(x["name"], x["tag"], float(d)), axis=1))
            f_ecif.close()
        print("\n- Finished -\n")

    def collect_ld(self, tag):
        """
        collect ligand descriptor data from sdf files
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        :return:
        """
        def calc_ld(iname, itag):
            if itag == "core":
                return
            if itag == "general":
                itag = "general-minus-refined"
            ligand_file = os.path.join(dataset_2016_dir[itag], iname, "{}_ligand.sdf".format(iname))
            ret = self._ecif_helper.get_ligand_features_by_file_sdf(ligand_file)
            if ret is not None:
                ld_writer.writerow([iname] + list(ret))
            return

        print("\nCollecting Ligand Descriptor data: \n")
        aff_data = pd.read_csv(os.path.join(data_dir, "affinity_data_2016.csv"))
        # for functional testing
        # aff_data = aff_data[0:50]
        if tag is None or tag == "":
            f_ld = open(os.path.join(meif_data_dir, "LD.csv"), 'a', newline='')
        else:
            f_ld = open(os.path.join(meif_data_dir, "LD_{}.csv".format(tag)), 'a', newline='')
        # write csv header
        ld_header = ["PDB"] + LIGAND_DESC
        ld_writer = csv.writer(f_ld)
        ld_writer.writerow(ld_header)
        list(aff_data.apply(lambda x: calc_ld(x["name"], x["tag"]), axis=1))
        f_ld.close()
        print("\n- Finished -\n")

    def joint_collect(self, tag):
        """
        jointly collect ecif, ld, and ground true value for each affinity entity
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        :return:
        """
        def calc_ecif(iname, itag, ipk, distance_cutoff):
            if itag == "general":
                itag = "general-minus-refined"

            if itag == "aug":
                sp = iname.split("_")
                target = sp[0]
                lig = sp[1]
                protein_file = os.path.join(data_aug_dir["proteins"], "{}.pdb".format(target))
                ligand_file = os.path.join(data_aug_dir["ligands"], target, "{}_ledock001.mol2".format(lig))
                try:
                    ecif_list = self._ecif_helper.get_ecif_mol2(protein_file,
                                                                ligand_file,
                                                                float(distance_cutoff))
                    ld_list = list(self._ecif_helper.get_ligand_features_by_file_mol2(ligand_file))
                except Exception:
                    return
            else:
                protein_file = os.path.join(dataset_2016_dir[itag], iname, "{}_protein.pdb".format(iname))
                ligand_file = os.path.join(dataset_2016_dir[itag], iname, "{}_ligand.sdf".format(iname))
                try:
                    ecif_list = self._ecif_helper.get_ecif(protein_file,
                                                           ligand_file,
                                                           float(distance_cutoff))
                    ld_list = list(self._ecif_helper.get_ligand_features_by_file_sdf(ligand_file))
                except Exception:
                    return

            if itag == "core":
                ecif_test_writer.writerow([iname] + ecif_list + ld_list + [ipk])
            else:
                ecif_train_writer.writerow([iname] + ecif_list + ld_list + [ipk])
            return

        print("\nCollecting ECIF data: \n")
        # affinity data are sorted by tag {core|general|refined}
        aff_data = pd.read_csv(aff_data_file).sort_values(by="tag", axis=0, ascending=True)
        # for functional testing
        # aff_data = aff_data[0:50]
        ecif_header = self._ecif_helper.get_possible_pl() + LIGAND_DESC
        fp_len = len(ecif_header)
        ecif_header = ["PDB"] + ecif_header + ["pk"]
        for d in self._distance_cutoffs:
            print("\n distance_cutoff: {}\n".format(d))
            if tag is None or tag == "":
                f_ecif_test = open(os.path.join(meif_data_dir, "ECIF_test_{}.csv".format(d)), 'a', newline='')
                f_ecif_train = open(os.path.join(meif_data_dir, "ECIF_train_{}.csv".format(d)), 'a', newline='')
            else:
                f_ecif_test = open(os.path.join(meif_data_dir, "ECIF_test_{}_{}.csv".format(d, tag)), 'a', newline='')
                f_ecif_train = open(os.path.join(meif_data_dir, "ECIF_train_{}_{}.csv".format(d, tag)), 'a', newline='')
            # write csv header
            ecif_test_writer = csv.writer(f_ecif_test)
            ecif_train_writer = csv.writer(f_ecif_train)
            ecif_test_writer.writerow(ecif_header)
            ecif_train_writer.writerow(ecif_header)
            aff_data.apply(lambda x: calc_ecif(x["name"], x["tag"], x["pk"], float(d)), axis=1)
            f_ecif_test.close()
            f_ecif_train.close()
        print("\n- Finished -\n")

    def joint_collect_wold(self, tag='wold'):
        """
        jointly collect ecif(without ld) and ground true value for each affinity entity
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        :return:
        """
        def calc_ecif(iname, itag, ipk, distance_cutoff):
            if itag == "general":
                itag = "general-minus-refined"
            protein_file = os.path.join(dataset_2016_dir[itag], iname, "{}_protein.pdb".format(iname))
            ligand_file = os.path.join(dataset_2016_dir[itag], iname, "{}_ligand.sdf".format(iname))

            try:
                ecif_list = self._ecif_helper.get_ecif(protein_file,
                                                       ligand_file,
                                                       float(distance_cutoff))
            except Exception:
                return

            if itag == "core":
                ecif_test_writer.writerow([iname] + ecif_list + [ipk])
            else:
                ecif_train_writer.writerow([iname] + ecif_list + [ipk])
            return

        print("\nCollecting ECIF data: \n")
        # affinity data are sorted by tag {core|general|refined}
        aff_data = pd.read_csv(os.path.join(data_dir, "affinity_data_2016.csv"))\
            .sort_values(by="tag", axis=0, ascending=True)
        # for functional testing
        # aff_data = aff_data[0:50]
        ecif_header = self._ecif_helper.get_possible_pl()
        ecif_header = ["PDB"] + ecif_header + ["pk"]
        for d in self._distance_cutoffs:
            print("\n distance_cutoff: {}\n".format(d))
            if tag is None or tag == "":
                f_ecif_test = open(os.path.join(meif_data_dir, "ECIF_test_{}.csv".format(d)), 'a', newline='')
                f_ecif_train = open(os.path.join(meif_data_dir, "ECIF_train_{}.csv".format(d)), 'a', newline='')
            else:
                f_ecif_test = open(os.path.join(meif_data_dir, "ECIF_test_{}_{}.csv".format(d, tag)), 'a', newline='')
                f_ecif_train = open(os.path.join(meif_data_dir, "ECIF_train_{}_{}.csv".format(d, tag)), 'a', newline='')
            # write csv header
            ecif_test_writer = csv.writer(f_ecif_test)
            ecif_train_writer = csv.writer(f_ecif_train)
            ecif_test_writer.writerow(ecif_header)
            ecif_train_writer.writerow(ecif_header)
            aff_data.apply(lambda x: calc_ecif(x["name"], x["tag"], x["pk"], float(d)), axis=1)
            f_ecif_test.close()
            f_ecif_train.close()
        print("\n- Finished -\n")


if __name__ == '__main__':
    ecif_collector = ECIF_Collector([6.0])
    ecif_collector.joint_collect("AUG")
