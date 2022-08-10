"""
platform: macOS
env: any
name: config.py
project configurations
"""
import os

configs_dir = os.path.abspath(os.path.dirname(__file__))
project_dir = os.path.split(configs_dir)[0]
tmp_dir = os.path.join(project_dir, "tmp")
data_dir = os.path.join(project_dir, "data")

if __name__ == '__main__':
    pass
