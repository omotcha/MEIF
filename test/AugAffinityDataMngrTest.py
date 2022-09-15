"""
platform: any
env: any
name: AugAffinityDataMngrTest.py
AugAffinityDataMngr Tester
"""
import os.path

from configs.config import *
from data.AugAffinityDataMngr import AugAffinityDataMngr
from util.ECIF import ECIF
import pandas as pd
from rdkit import Chem

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


def test():
    f_prot = os.path.join(data_aug_dir["proteins"], "966c.pdb")
    f_ligd = os.path.join(data_aug_dir["ligands"], "966c", "83_ledock001.mol2")
    # ecif_helper = ECIF()
    # ecif_list = ecif_helper.get_ecif_mol2(f_prot, f_ligd, float(6.0))
    # ld_list = ecif_helper.get_ligand_features_by_file_mol2(f_ligd)
    m = Chem.MolFromMol2File(f_ligd, sanitize=False)
    print(m)
    # print(ecif_list)
    # print(ld_list)


if __name__ == '__main__':
    test()
