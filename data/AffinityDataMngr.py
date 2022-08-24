"""
platform: any
env: any
name: AffinityDataMngr.py
create affinity data file from dataset
"""
import os
from configs.config import dataset_version, platform_delimiter, index_file_dir, data_dir
import re
import pandas as pd


class AffinityDataMngr:
    _ds_version = 2020

    def __init__(self):
        self._ds_version = dataset_version

    # callables:
    def set_ds_version(self, ds_version):
        """
        pdbbind dataset version setter
        :param ds_version:
        :return:
        """
        if ds_version in (2013, 2016, 2019, 2020):
            self._ds_version = ds_version
        else:
            print("Error: dataset version can only be chosen from 2013|2016|2019|2020")

    def get_aff_data_from_index_file(self, f_index, tag):
        """
        store affinity data from index file to a dataframe
        :param f_index: the index file, something like INDEX_XXXXX_data.YYYY
        :param tag: tag to identify affinity data, if not set, it would be filled with index file name XXXXX
        :return:
        """
        names = []
        pks = []
        if tag is None:
            tag = f_index.split(platform_delimiter)[-1][6:-5].split("_")[0]
        with open(f_index, "r") as f:
            lines = f.readlines()
            for line in lines:
                seg = re.split(' +', line)
                if seg[0] != '#':
                    names.append(seg[0])
                    pks.append(seg[3])
        tags = [tag] * len(names)
        aff = pd.DataFrame({"name": names, "pk": pks, "tag": tags}, columns=["name", "pk", "tag"])
        return aff

    def write_aff_data(self):
        """
        save all affinity data to affinity_data.csv
        :return:
        """
        general_aff = self.get_aff_data_from_index_file(index_file_dir["general"], "general")
        refined_aff = self.get_aff_data_from_index_file(index_file_dir["refined"], "refined")
        core_aff = self.get_aff_data_from_index_file(index_file_dir["core"], "core")
        mask_gr = general_aff["name"].isin(refined_aff["name"])
        general_aff.loc[mask_gr, "tag"] = "refined"
        mask_gc = general_aff["name"].isin(core_aff["name"])
        general_aff.loc[mask_gc, "tag"] = "core"
        general_aff.to_csv(os.path.join(data_dir, "affinity_data.csv"))

    def get_aff_data(self):
        """
        get all affinity data in pandas dataframe format
        :return:
        """
        general_aff = self.get_aff_data_from_index_file(index_file_dir["general"], "general")
        refined_aff = self.get_aff_data_from_index_file(index_file_dir["refined"], "refined")
        core_aff = self.get_aff_data_from_index_file(index_file_dir["core"], "core")
        mask_gr = general_aff["name"].isin(refined_aff["name"])
        general_aff.loc[mask_gr, "tag"] = "refined"
        mask_gc = general_aff["name"].isin(core_aff["name"])
        general_aff.loc[mask_gc, "tag"] = "core"
        return general_aff
