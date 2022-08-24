"""
platform: any
env: any
name: meif_collector.py
collect MEIF fingerprints of all protein-ligand pairs
"""
import os.path

import pandas as pd

from configs.config import data_dir, meif_data_dir, dataset_dir
from util.MEIF import MEIF, LIGAND_DESC
import csv


class MEIF_Collector:
    _distance_cutoffs = []
    _meif_helper = None

    def __init__(self, dc_list):
        if dc_list is None or len(dc_list) == 0:
            self._distance_cutoffs = [6.0]
        else:
            self._distance_cutoffs = dc_list
        self._meif_helper = MEIF()

    # callables

    def collect_meif(self, tag):
        """
        collect MEIF data from pdb files and sdf files
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        :return:
        """
        def calc_meif(iname, itag, distance_cutoff):
            if itag == "general":
                itag = "general-minus-refined"
            protein_file = os.path.join(dataset_dir[itag], iname, "{}_protein.pdb".format(iname))
            ligand_file = os.path.join(dataset_dir[itag], iname, "{}_ligand.sdf".format(iname))
            # suppose all cutoffs are changing in same order, this pattern can be modified anyway
            meif_list = self._meif_helper.get_meif(protein_file,
                                                   ligand_file,
                                                   float(distance_cutoff),
                                                   float(distance_cutoff),
                                                   float(distance_cutoff),
                                                   None)
            meif_writer.writerow([iname] + meif_list)
            return

        print("\nCollecting MEIF data: \n")
        aff_data = pd.read_csv(os.path.join(data_dir, "affinity_data.csv"))
        # for functional testing
        # aff_data = aff_data[0:50]
        for d in self._distance_cutoffs:
            print("\n distance_cutoff: {}\n".format(d))
            if tag is None or tag == "":
                f_meif = open(os.path.join(meif_data_dir, "MEIF_{}.csv".format(d)), 'a', newline='')
            else:
                f_meif = open(os.path.join(meif_data_dir, "MEIF_{}_{}.csv".format(d, tag)), 'a', newline='')
            # write csv header
            meif_header = ["PDB"] + self._meif_helper.get_elem_head()
            meif_writer = csv.writer(f_meif)
            meif_writer.writerow(meif_header)
            list(aff_data.apply(lambda x: calc_meif(x["name"], x["tag"], float(d)), axis=1))
            f_meif.close()
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
            ld_list = list(self._meif_helper.get_ligand_features_by_file(ligand_file))
            ld_writer.writerow([iname] + ld_list)
            return

        print("\nCollecting Ligand Descriptor data: \n")
        aff_data = pd.read_csv(os.path.join(data_dir, "affinity_data.csv"))
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


if __name__ == '__main__':
    meif_collector = MEIF_Collector([6.0])
    meif_collector.collect_meif("")
