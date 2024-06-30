import os
import shutil
from pathlib import Path

DEFAULT_DATA_FOLDER = str(Path(__file__).parents[1]) + '/data'
RBC2_DATA_PATH = os.getenv('RBC2_DATA_PATH')

if RBC2_DATA_PATH is None:
    # if environment variable not set, use default path
    path_to_data_folder = DEFAULT_DATA_FOLDER

else:
    path_to_data_folder = RBC2_DATA_PATH

    # copy the existing data folder to the env path if it doesn't exist already
    if not os.path.exists(path_to_data_folder):
        shutil.copytree(DEFAULT_DATA_FOLDER, path_to_data_folder)


if __name__ == '__main__':
    print(path_to_data_folder)
    print(DEFAULT_DATA_FOLDER)
    print(RBC2_DATA_PATH)



