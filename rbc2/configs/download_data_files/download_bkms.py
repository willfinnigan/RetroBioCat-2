import os
from pathlib import Path

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_file, download_gdrive_file


def does_bkms_exist() -> bool:
    if not os.path.exists(f"{path_to_data_folder}/bkms/policy_model/variables/variables.data-00000-of-00001"):
        return False

    if not os.path.exists(f"{path_to_data_folder}/bkms/bkms_metadata.hdf"):
        return False

    if not os.path.exists(f"{path_to_data_folder}/bkms/bkms_templates_only.hdf"):
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

    directory = f"{path_to_data_folder}/bkms/"
    file_id = "1Cgufom9YA3ElYO7v9jy6Vi2qWiWOMvHX"
    filename = "bkms_metadata.hdf"
    filepath = f"{directory}/{filename}"
    download_gdrive_file(file_id, filepath)

    file_id = "1CdHaaqfWzM70e-csxiZA6cdWEwpY9HA5"
    filename = "bkms_templates_only.hdf"
    filepath = f"{directory}/{filename}"
    download_gdrive_file(file_id, filepath)





if __name__ == '__main__':
    download_bkms_model()