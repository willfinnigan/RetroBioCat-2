import os
from pathlib import Path

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_file


def does_bkms_exist() -> bool:
    if not os.path.exists(f"{path_to_data_folder}/bkms/policy_model/variables/variables.data-00000-of-00001"):
        return False

    return True

def download_bkms_model():
    bkms_model = "https://github.com/itai-levin/bkms-data/blob/master/models/bkms/1/variables/variables.data-00000-of-00001?raw=true"

    # if folder doesn't exist, create it with Pathlib
    directory = f"{path_to_data_folder}/bkms/policy_model/variables"

    Path(directory).mkdir(parents=True, exist_ok=True)

    filename = "variables.data-00000-of-00001"
    filepath = f"{directory}/{filename}"
    download_file(bkms_model, filepath)



if __name__ == '__main__':
    download_bkms_model()