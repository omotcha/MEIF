"""
platform: any
env: any
name: MEIFCollectorTest.py
MEIF Collector Tester
"""
from configs.config import *
from preprocess.meif_collector import MEIF_Collector


def testCollectMEIF():
    meif_collector = MEIF_Collector([6.0])
    meif_collector.collect_meif("")


def testCollectLD():
    ld_collector = MEIF_Collector([6.0])
    ld_collector.collect_ld("")


if __name__ == '__main__':
    testCollectLD()
