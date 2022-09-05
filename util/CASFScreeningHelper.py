"""
platform: any
env: any
name: CASFScreeningHelper.py
CASF screening analysis supporter
"""
import os
from configs.config import casf_dir
import math
import multiprocessing
import time
import pandas as pd


def target_based_data_collector_worker(targets, d_input, d_output):
    """

    :param targets:
    :param d_input:
    :param d_output:
    :return:
    """
    for tid in targets:
        target_dfs = []
        for f in os.listdir(d_input):
            if f.startswith(tid):
                target_dfs.append(pd.read_csv(os.path.join(d_input, f)))
        df = pd.concat(target_dfs, ignore_index=True)
        df.to_csv(os.path.join(d_output, "{}_score.dat".format(tid)))



def target_based_data_collector_modulator(n_workers=4, d_input=None, d_output=None):
    """
    for CASF screening power analysis, this function serves as a data collector for VS results
    :param n_workers:
    :param d_input: input directory of data files "{target}_{ligand}_score.dat"
    :param d_output: output directory of data files "{target}_score.dat"
    :return:
    """
    def get_targets():
        ret = []
        with open(os.path.join(casf_dir["target_info"]), 'r') as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    ret.append(line.split(" ")[0])
        return ret
    targets = get_targets()
    size = math.ceil(len(targets) / n_workers)
    worker_tasks = [targets[i:i + size] for i in range(0, len(targets), size)]
    start = time.perf_counter()
    pool = multiprocessing.Pool(n_workers)
    pool.starmap(target_based_data_collector_worker, zip(worker_tasks, [d_input] * n_workers, [d_output] * n_workers))
    end = time.perf_counter()
    print('\n')
    print('run time: {} seconds'.format(round(end - start)))
