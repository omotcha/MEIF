"""
platform: any
env: any
name: RDKitHelperTest.py
RDKitHelper Tester
"""

import os
from util.RDKitHelper import get_decoy_names
from configs.config import casf_dir


def testGetDecoyNames():
    names = get_decoy_names(os.path.join(casf_dir["decoys"], "1a30_decoys.mol2"))
    print(names)


if __name__ == '__main__':
    testGetDecoyNames()
