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


if __name__ == '__main__':
    testGetPLPairWithDecoysCached()
