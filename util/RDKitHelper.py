"""
platform: win
env: any with RDKit
name: RDKitHelper.py
extend RDKit functionality
"""
import os
from rdkit import Chem


def Mol2MolSupplier(file=None, sanitize=True):
    """
    read mol2 file with multiple molecules
    modified from: https://chem-workflows.com/articles/2019/07/18/building-a-multi-molecule-mol2-reader-for-rdkit/
    :param file: mol2 file
    :param sanitize:
    :return:
    """
    mols = []
    with open(file, 'r') as f:
        line = f.readline()
        while not line.startswith("@<TRIPOS>MOLECULE"):
            line = f.readline()
        m = None
        while not f.tell() == os.fstat(f.fileno()).st_size:
            if line.startswith("@<TRIPOS>MOLECULE"):
                mol = [line]
                line = f.readline()
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()
                    if f.tell() == os.fstat(f.fileno()).st_size:
                        mol.append(line)
                        break
                mol[-1] = mol[-1].rstrip()  # removes blank line at file end
                block = ",".join(mol).replace(',', '')
                m = Chem.MolFromMol2Block(block, sanitize=sanitize)
            mols.append(m)
    return mols
