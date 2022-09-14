"""
platform: any
env: any
name: AugAffinityDataMngrTest.py
AugAffinityDataMngr Tester
"""

from configs.config import *
from data.AugAffinityDataMngr import AugAffinityDataMngr
import pandas as pd

affdata_helper = AugAffinityDataMngr()


def testFather():
    print("\n")
    core_aff = affdata_helper.get_aff_data_from_index_file(index_file_2016_dir["core"], "core")
    print(core_aff)


def testSon():
    print("\n")
    aug_aff = affdata_helper.get_aff_data_from_augmented_label_file(data_aug_dir["label"])
    print(aug_aff)


def testWriteAll():
    affdata_helper.write_aug_aff_data()


if __name__ == '__main__':
    testWriteAll()
