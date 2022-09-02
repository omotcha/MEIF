"""
platform: any
env: any
name: ecifp_collector.py
collect ECIFP fingerprints of all protein-ligand pairs
"""
import os.path

import pandas as pd
from configs.config import data_dir, meif_data_dir, dataset_dir
from util.ECIFP import ECIFP, LIGAND_DESC
import csv


class ECIFP_Collector:
    _distance_cutoffs = []
    _ecifp_helper = None

    def __init__(self, dc_list):
        if dc_list is None or len(dc_list) == 0:
            self._distance_cutoffs = [6.0]
        else:
            self._distance_cutoffs = dc_list
        self._ecifp_helper = ECIFP()

    # callables

    def collect_ecifp(self, tag):
        """
        collect ECIFP data from pdb files and sdf files
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        :return:
        """
        def calc_ecifp(iname, itag, distance_cutoff):

            if itag == "general":
                itag = "general-minus-refined"
            protein_file = os.path.join(dataset_dir[itag], iname, "{}_protein.pdb".format(iname))
            ligand_file = os.path.join(dataset_dir[itag], iname, "{}_ligand.sdf".format(iname))
            ecifp_list = self._ecifp_helper.get_ecifp(protein_file,
                                                      ligand_file,
                                                      float(distance_cutoff))
            if itag == "core":
                ecifp_test_writer.writerow([iname] + ecifp_list)
            else:
                ecifp_train_writer.writerow([iname] + ecifp_list)
            return

        print("\nCollecting ECIFP data: \n")
        aff_data = pd.read_csv(os.path.join(data_dir, "affinity_data.csv"))
        # for functional testing
        # aff_data = aff_data[0:50]
        for d in self._distance_cutoffs:
            print("\n distance_cutoff: {}\n".format(d))
            if tag is None or tag == "":
                f_ecifp_test = open(os.path.join(meif_data_dir, "ECIFP_test_{}.csv".format(d)), 'a', newline='')
                f_ecifp_train = open(os.path.join(meif_data_dir, "ECIFP_train_{}.csv".format(d)), 'a', newline='')
            else:
                f_ecifp_test = open(os.path.join(meif_data_dir, "ECIFP_test_{}_{}.csv".format(d, tag)), 'a', newline='')
                f_ecifp_train = open(os.path.join(meif_data_dir, "ECIFP_train_{}_{}.csv".format(d, tag)), 'a', newline='')
            # write csv header
            ecifp_header = ["PDB"] + self._ecifp_helper.get_possible_pl()
            ecifp_test_writer = csv.writer(f_ecifp_test)
            ecifp_train_writer = csv.writer(f_ecifp_train)
            ecifp_test_writer.writerow(ecifp_header)
            ecifp_train_writer.writerow(ecifp_header)
            list(aff_data.apply(lambda x: calc_ecifp(x["name"], x["tag"], float(d)), axis=1))
            f_ecifp_test.close()
            f_ecifp_train.close()
        print("\n- Finished -\n")

    def collect_ld(self, tag):
        """
        collect ligand descriptor data from sdf files
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        :return:
        """
        def calc_ld(iname, itag):
            if itag == "general":
                itag = "general-minus-refined"
            ligand_file = os.path.join(dataset_dir[itag], iname, "{}_ligand.sdf".format(iname))
            ret = self._ecifp_helper.get_ligand_features_by_file(ligand_file)
            if ret is not None:
                if itag == "core":
                    ld_test_writer.writerow([iname] + list(ret))
                else:
                    ld_train_writer.writerow([iname] + list(ret))
            return

        print("\nCollecting Ligand Descriptor data: \n")
        aff_data = pd.read_csv(os.path.join(data_dir, "affinity_data.csv"))
        # for functional testing
        # aff_data = aff_data[0:50]
        if tag is None or tag == "":
            f_ld_test = open(os.path.join(meif_data_dir, "LD_test.csv"), 'a', newline='')
            f_ld_train = open(os.path.join(meif_data_dir, "LD_train.csv"), 'a', newline='')
        else:
            f_ld_test = open(os.path.join(meif_data_dir, "LD_test_{}.csv".format(tag)), 'a', newline='')
            f_ld_train = open(os.path.join(meif_data_dir, "LD_train_{}.csv".format(tag)), 'a', newline='')
        # write csv header
        ld_header = ["PDB"] + LIGAND_DESC
        ld_test_writer = csv.writer(f_ld_test)
        ld_train_writer = csv.writer(f_ld_train)
        ld_test_writer.writerow(ld_header)
        ld_train_writer.writerow(ld_header)
        list(aff_data.apply(lambda x: calc_ld(x["name"], x["tag"]), axis=1))
        f_ld_test.close()
        f_ld_train.close()
        print("\n- Finished -\n")

    def joint_collect(self, tag):
        """
        jointly collect ecifp, ld, and ground true value for each affinity entity
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        :return:
        """
        def calc_ecifp(iname, itag, ipk, distance_cutoff):
            if itag == "general":
                itag = "general-minus-refined"
            protein_file = os.path.join(dataset_dir[itag], iname, "{}_protein.pdb".format(iname))
            ligand_file = os.path.join(dataset_dir[itag], iname, "{}_ligand.sdf".format(iname))

            try:
                ecifp_list = self._ecifp_helper.get_ecifp(protein_file,
                                                          ligand_file,
                                                          float(distance_cutoff))
                ld_list = list(self._ecifp_helper.get_ligand_features_by_file(ligand_file))
            except Exception:
                return

            if itag == "core":
                ecifp_test_writer.writerow([iname] + ecifp_list + ld_list + [ipk])
            else:
                ecifp_train_writer.writerow([iname] + ecifp_list + ld_list + [ipk])
            return

        print("\nCollecting ECIFP data: \n")
        # affinity data are sorted by tag {core|general|refined}
        aff_data = pd.read_csv(os.path.join(data_dir, "affinity_data.csv")).sort_values(by="tag", axis=0, ascending=True)
        # for functional testing
        # aff_data = aff_data[0:50]
        ecifp_header = self._ecifp_helper.get_possible_pl() + LIGAND_DESC
        ecifp_header = ["PDB"] + ecifp_header + ["pk"]
        for d in self._distance_cutoffs:
            print("\n distance_cutoff: {}\n".format(d))
            if tag is None or tag == "":
                f_ecifp_test = open(os.path.join(meif_data_dir, "ECIFP_test_{}.csv".format(d)), 'a', newline='')
                f_ecifp_train = open(os.path.join(meif_data_dir, "ECIFP_train_{}.csv".format(d)), 'a', newline='')
            else:
                f_ecifp_test = open(os.path.join(meif_data_dir, "ECIFP_test_{}_{}.csv".format(d, tag)), 'a', newline='')
                f_ecifp_train = open(os.path.join(meif_data_dir, "ECIFP_train_{}_{}.csv".format(d, tag)), 'a', newline='')
            # write csv header
            ecifp_test_writer = csv.writer(f_ecifp_test)
            ecifp_train_writer = csv.writer(f_ecifp_train)
            ecifp_test_writer.writerow(ecifp_header)
            ecifp_train_writer.writerow(ecifp_header)
            aff_data.apply(lambda x: calc_ecifp(x["name"], x["tag"], x["pk"], float(d)), axis=1)
            f_ecifp_test.close()
            f_ecifp_train.close()
        print("\n- Finished -\n")


if __name__ == '__main__':
    ecifp_collector = ECIFP_Collector([6.0])
    ecifp_collector.joint_collect("")
