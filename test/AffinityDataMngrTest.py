"""
platform: any
env: any
name: AffinityDataMngrTest.py
AffinityDataMngr Tester
"""

from configs.config import *
from data.AffinityDataMngr import AffinityDataMngr
import pandas as pd

affdata_helper = AffinityDataMngr()


def testDataSubsetRelation():
    print("\n")

    general_aff = affdata_helper.get_aff_data_from_index_file(index_file_dir["general"], "general")
    refined_aff = affdata_helper.get_aff_data_from_index_file(index_file_dir["refined"], "refined")
    core_aff = affdata_helper.get_aff_data_from_index_file(index_file_dir["core"], "core")
    core_in_general = core_aff[core_aff["name"].isin(general_aff["name"])]
    core_in_refined = core_aff[core_aff["name"].isin(refined_aff["name"])]
    refined_in_general = refined_aff[refined_aff["name"].isin(general_aff["name"])]
    assert core_in_general.shape == core_aff.shape
    assert core_in_refined.shape <= core_aff.shape
    assert refined_in_general.shape == refined_aff.shape


def testAllAffData():
    print("\n")
    affdata_helper.write_aff_data()


if __name__ == '__main__':
    testAllAffData()
