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
            if m is not None:
                mols.append(m)
    return mols


def Mol2MolSupplier_test(file=None, sanitize=True):
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
                if m is not None:
                    mols.append(m)
                else:
                    print("\n")
                    print(block)
    return mols


def get_decoy_names_mol2(file=None):
    """
    get the decoy names from mol2 file with multiple decoys
    :param file: mol2 file
    :return:
    """
    names = []
    with open(file, 'r') as f:
        line = f.readline()
        while not line.startswith("@<TRIPOS>MOLECULE"):
            line = f.readline()
        while not f.tell() == os.fstat(f.fileno()).st_size:
            if line.startswith("@<TRIPOS>MOLECULE"):
                line = f.readline()
                names.append(line[:-1])
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    line = f.readline()
                    if f.tell() == os.fstat(f.fileno()).st_size:
                        break

    return names


def get_decoy_names_sdf(file=None):
    """
    get the decoy names from sdf file with multiple decoys
    :param file: sdf file
    :return:
    """
    names = []
    decoys = Chem.SDMolSupplier(file)
    for decoy in decoys:
        names.append(decoy.GetProp('_Name'))
    return names
