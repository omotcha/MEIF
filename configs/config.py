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

# dataset used
dataset_version = 2020
dataset_dir = {
    "base": "E:\\datasets\\dataset\\pdbbind2020",
    "general-minus-refined": "E:\\datasets\\dataset\\pdbbind2020\\general-minus-refined",
    "refined": "E:\\datasets\\dataset\\pdbbind2020\\refined-set",
    "core": "E:\\datasets\\dataset\\pdbbind2016\\coreset"
}
index_file_dir = {
    "general": "E:\\datasets\\dataset\\pdbbind2020\\index\\INDEX_general_PL_data.2020",
    "refined": "E:\\datasets\\dataset\\pdbbind2020\\index\\INDEX_refined_data.2020",
    "core": "E:\\datasets\\dataset\\pdbbind2016\\PDBbind_2016_plain_text_index\\index\\INDEX_core_data.2016"
}


if __name__ == '__main__':
    print(project_dir)
