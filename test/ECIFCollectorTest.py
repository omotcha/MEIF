"""
platform: any
env: any
name: ECIFCollectorTest.py
ECIF Collector Tester
"""
from configs.config import *
from util.ECIF import ECIF


def testMol2Collection():
    f_pro = os.path.join(data_aug_dir["proteins"], "1cr6.pdb")
    f_lig = os.path.join(data_aug_dir["ligands"], "1cr6", "1_ledock001.mol2")
    ecif_helper = ECIF()
    ecif = ecif_helper.get_ecif_mol2(f_pro, f_lig, float(6.0))
    ld = ecif_helper.get_ligand_features_by_file_mol2(f_lig)
    print(ecif + list(ld))


if __name__ == '__main__':
    testMol2Collection()
