"""
platform: win
env: any
name: casf_support.py
CASF analysis supporter
"""
import os
from configs.config import casf_dir
from util.RDKitHelper import Mol2MolSupplier
import pandas as pd
from rdkit.Chem import AllChem


def get_decoys(lid):
    """
    get all the decoys from given ligand id
    :param lid: ligand id
    :return:
    """
    f_decoy = os.path.join(casf_dir["decoys"], "{}_decoys.mol2".format(lid))
    # f_single_lig_test = os.path.join(casf_dir["core"], lid, "{}_ligand.mol2".format(lid))
    # ligands = Chem.MolFromMol2File(f_single_lig_test)
    ligands = Mol2MolSupplier(f_decoy, sanitize=True)
    print(len(ligands))
    print(ligands[0].GetNumAtoms())


if __name__ == '__main__':
    get_decoys("1a30")
