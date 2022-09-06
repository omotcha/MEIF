"""
platform: any
env: any
name: RDKitHelperTest.py
RDKitHelper Tester
"""

import os
from util.RDKitHelper import Mol2MolSupplier, get_decoy_names_mol2
from rdkit import Chem
from configs.config import casf_dir, tmp_dir


def testGetDecoyNames():
    names = get_decoy_names_mol2(os.path.join(casf_dir["decoys"], "1bzc_decoys.mol2"))
    print(names)


def testGetDecoys():
    decoys = Mol2MolSupplier(os.path.join(casf_dir["decoys"], "1bzc_decoys.mol2"), sanitize=True)
    print(len(decoys))


def testGetDecoysFromSDF(f):
    decoys = Chem.SDMolSupplier(f)
    for decoy in decoys:
        print(decoy.GetProp('_Name'))
    print(len(decoys))


if __name__ == '__main__':
    testGetDecoysFromSDF(os.path.join(tmp_dir, "1c5z_decoys.sdf"))
