"""
platform: any
env: any
name: ECIFPCollectorTest.py
ECIFP Collector Tester
"""
from configs.config import *
from preprocess.ecifp_collector import ECIFP_Collector


def testCollectECIFP():
    ecifp_collector = ECIFP_Collector([6.0])
    ecifp_collector.joint_collect("")


# def testCollectLD():
#     ld_collector = ECIFP_Collector([6.0])
#     ld_collector.collect_ld("")


if __name__ == '__main__':
    testCollectECIFP()
