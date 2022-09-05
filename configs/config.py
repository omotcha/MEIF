"""
platform: win
env: any
name: config.py
project configurations
"""
import os

# platform specific
platform_delimiter = "\\"

# project structure
configs_dir = os.path.abspath(os.path.dirname(__file__))
project_dir = os.path.split(configs_dir)[0]
tmp_dir = os.path.join(project_dir, "tmp")
data_dir = os.path.join(project_dir, "data")
meif_data_dir = os.path.join(project_dir, "meif_data")
test_dir = os.path.join(project_dir, "test")
model_test_dir = os.path.join(project_dir, "model_test")
log_dir = os.path.join(project_dir, "log")

# dataset used
dataset_version = 2020
dataset_dir = {
    "base": "E:\\datasets\\dataset\\pdbbind2020",
    "general-minus-refined": "E:\\datasets\\dataset\\pdbbind2020\\general-minus-refined",
    "refined": "E:\\datasets\\dataset\\pdbbind2020\\refined-set",
    "core": "E:\\datasets\\dataset\\pdbbind2016\\coreset"
}

dataset_2016_dir = {
    "base": "E:\\datasets\\dataset\\pdbbind2016",
    "general-minus-refined": "E:\\datasets\\dataset\\pdbbind2016\\general-set-except-refined",
    "refined": "E:\\datasets\\dataset\\pdbbind2020\\refined-set",
    "core": "E:\\datasets\\dataset\\pdbbind2016\\coreset"
}

index_file_dir = {
    "general": "E:\\datasets\\dataset\\pdbbind2020\\index\\INDEX_general_PL_data.2020",
    "refined": "E:\\datasets\\dataset\\pdbbind2020\\index\\INDEX_refined_data.2020",
    "core": "E:\\datasets\\dataset\\pdbbind2016\\PDBbind_2016_plain_text_index\\index\\INDEX_core_data.2016"
}

index_file_2016_dir = {
    "general": "E:\\datasets\\dataset\\pdbbind2016\\PDBbind_2016_plain_text_index\\index\\INDEX_general_PL_data.2016",
    "refined": "E:\\datasets\\dataset\\pdbbind2016\\PDBbind_2016_plain_text_index\\index\\INDEX_refined_data.2016",
    "core": "E:\\datasets\\dataset\\pdbbind2016\\PDBbind_2016_plain_text_index\\index\\INDEX_core_data.2016"
}

casf_dir = {
    "base": "E:\\datasets\\casf2016\\CASF-2016",
    "core": "E:\\datasets\\casf2016\\CASF-2016\\coreset",
    "decoys": "E:\\datasets\\casf2016\\CASF-2016\\decoys_docking",
    "decoys_sdf": "E:\\datasets\\casf2016\\CASF-2016\\decoys_docking_sdf",
    "target_info": "E:\\datasets\\casf2016\\CASF-2016\\power_screening\\TargetInfo.dat",
    "decoys_screening": "E:\\datasets\\casf2016\\CASF-2016\\decoys_screening",
    "decoys_screening_sdf": "E:\\datasets\\casf2016\\CASF-2016\\decoys_screening_sdf"
}

# models
ecif_base = data_dir

ecif_gbt = os.path.join(ecif_base, "ECIF6_LD_GBT.pkl")
ecif_ag_base = os.path.join(project_dir, "train", "AutogluonModels", "ag-20220902_034242", "models")
ecif_catboost = os.path.join(ecif_ag_base, "CatBoost", "model.pkl")
ecif_lightgbmxt = os.path.join(ecif_ag_base, "LightGBMXT", "model.pkl")

ecifp_base = os.path.join(project_dir, "train", "AutogluonModels", "ag-20220902_011835", "models")
ecifp_gbt = os.path.join(ecif_base, "ECIFP_LD_GBT.pkl")
ecifp_catboost = os.path.join(ecifp_base, "CatBoost", "model.pkl")
ecifp_lightgbm = os.path.join(ecifp_base, "LightGBM", "model.pkl")


if __name__ == '__main__':
    print(project_dir)
