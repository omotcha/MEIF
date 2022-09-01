"""
platform: any
env: any
name: ECIFPDockingPowerTest.py
ECIFP Docking Power Tester
"""
import os
from configs.config import casf_dir
from util.ECIFP import ECIFP


def testLoadLigands():
    print("\n")
    helper = ECIFP()
    lid = "1a30"
    f_decoy = os.path.join(casf_dir["decoys"], "{}_decoys.mol2".format(lid))
    dfs = helper.load_ligands(f_decoy)
    print(len(dfs))
    for df in dfs[:10]:
        print(df)


def testGetPLPairWithDecoysCached():
    print("\n")
    helper = ECIFP()
    lid = "1a30"
    f_prot = os.path.join(casf_dir["core"], lid, "{}_protein.pdb".format(lid))
    f_decoy = os.path.join(casf_dir["decoys"], "{}_decoys.mol2".format(lid))

    helper.cache_protein(f_prot)
    pairs = helper.get_pl_pairs_with_decoys_cached(f_decoy, 6.0)
    print(len(pairs))
    for df in pairs[:10]:
        print(df)


def testDecoysLD():
    print("\n")
    helper = ECIFP()
    lid = "1a30"
    f_decoy = os.path.join(casf_dir["decoys"], "{}_decoys.mol2".format(lid))
    lds = helper.test_ligand_features_decoys(f_decoy)
    same = True
    # for i in range(len(lds)-1):
    #     if lds[0] != lds[i+1]:
    #         same = False
    #         print("differs on ID: {}".format(i+1))
    #
    # print(same)

    for i in range(len(lds[0])):
        if lds[0][i] != lds[30][i]:
            print(i)

    print(lds[0])
    print(lds[30])


if __name__ == '__main__':
    testDecoysLD()
