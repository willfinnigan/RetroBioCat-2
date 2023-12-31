import os
import shutil
from pathlib import Path

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_file


def does_aizynthfinder_exist() -> bool:
    if not os.path.exists(f"{path_to_data_folder}/aizynthfinder/uspto_model.hdf5"):
        return False
    if not os.path.exists(f"{path_to_data_folder}/aizynthfinder/uspto_templates.hdf5"):
        return False
    return True

def download_aizynthfinder_model():
    aizynthfinder_model = "https://figshare.com/ndownloader/files/23086454"
    aizynthfinder_templates = "https://figshare.com/ndownloader/files/23086457"

    # if aizynthfinder folder doesn't exist, create it with Pathlib
    directory = f"{path_to_data_folder}/aizynthfinder"

    Path(directory).mkdir(parents=True, exist_ok=True)

    filename = "uspto_model.hdf5"
    filepath = f"{directory}/{filename}"
    download_file(aizynthfinder_model, filepath)

    filename = "uspto_templates.hdf5"
    filepath = f"{directory}/{filename}"
    download_file(aizynthfinder_templates, filepath)



if __name__ == '__main__':
    download_aizynthfinder_model()

