import os
import pandas as pd

final_kikd_result = "/dssg/home/acct-ins-hl/ins-hl/chaohao/Database/BindingDB/KiKd/KiKd_final_uniqPDBID.tsv"
docking_set = "/dssg/home/acct-ins-hl/ins-hl/hlcs/hjw/AutoLeDock/data/docking_set"

# final_kikd_result = "E:\\datasets\\data_aug\\KiKd_final_uniqPDBID.tsv"
# docking_set = "E:\\datasets\\data_aug\\docking_set"


def labeling():
    kikd_data = pd.read_csv(final_kikd_result, delimiter='\t', header=None)\
        .rename(columns={0: "raw_id", 1: "pk", 2: "smiles", 3: "target"})
    # kikd_data["docking_id"] = 0
    all_results = pd.DataFrame(columns=["docking_id", "smiles", "docking_score", "target"])
    for target in os.listdir(docking_set):
        result = os.path.join(docking_set, target, "task_1", "task_1_result.csv")
        if not os.path.isfile(result):
            continue
        result_data = pd.read_csv(result, header=None)\
            .rename(columns={0: "docking_id", 1: "smiles", 2: "docking_score"})
        result_data["target"] = target.upper()
        all_results = all_results.append(result_data)
    kikd_data = kikd_data.merge(all_results, left_on=["smiles", "target"], right_on=["smiles", "target"])\
        .drop(columns=["raw_id", "smiles", "docking_score"])
    kikd_data.to_csv(os.path.join(docking_set, "labeled_kikd_data.csv"),index=False)


if __name__ == '__main__':
    labeling()
