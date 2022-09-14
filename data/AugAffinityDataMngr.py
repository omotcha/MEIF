"""
platform: any
env: any
name: AugAffinityDataMngr.py
create augmented affinity data file from augmented dataset
"""

import os
from configs.config import index_file_dir, index_file_2016_dir, data_aug_dir, data_dir
from data.AffinityDataMngr import AffinityDataMngr
import pandas as pd


class AugAffinityDataMngr(AffinityDataMngr):

    def get_aff_data_from_augmented_label_file(self, aug_lb_f):
        """
        store affinity data from augmented label file to a dataframe
        :param aug_lb_f: augmented label file
        :return:
        """
        aug_data = pd.read_csv(aug_lb_f)
        aug_data["name"] = (aug_data["target"].astype(str).str.lower() + "_" + aug_data["docking_id"].astype(str))
        aug_data = aug_data.drop(columns=["target", "docking_id"])
        aug_data["tag"] = "aug"
        return aug_data

    def write_aug_aff_data(self):
        if self._ds_version == 2016:
            general_aff = self.get_aff_data_from_index_file(index_file_2016_dir["general"], "general")
            refined_aff = self.get_aff_data_from_index_file(index_file_2016_dir["refined"], "refined")
            core_aff = self.get_aff_data_from_index_file(index_file_2016_dir["core"], "core")
        else:
            general_aff = self.get_aff_data_from_index_file(index_file_dir["general"], "general")
            refined_aff = self.get_aff_data_from_index_file(index_file_dir["refined"], "refined")
            core_aff = self.get_aff_data_from_index_file(index_file_dir["core"], "core")
        mask_gr = general_aff["name"].isin(refined_aff["name"])
        general_aff.loc[mask_gr, "tag"] = "refined"
        mask_gc = general_aff["name"].isin(core_aff["name"])
        general_aff.loc[mask_gc, "tag"] = "core"

        general_aff = general_aff.append(self.get_aff_data_from_augmented_label_file(data_aug_dir["label"]))

        if self._ds_version == 2016:
            general_aff.to_csv(os.path.join(data_dir, "aug_affinity_data_2016.csv"), index=False)
        else:
            general_aff.to_csv(os.path.join(data_dir, "aug_affinity_data.csv"), index=False)


