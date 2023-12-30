import os
from pathlib import Path

# get path_to_data_folder frm enviorment variable
path_to_data_folder = os.getenv('RBC2_DATA_PATH')

if path_to_data_folder is None:
    # if enviorment variable not set, use default path
    path_to_data_folder = str(Path(__file__).parents[1]) + '/data'

if __name__ == '__main__':
    print(path_to_data_folder)

