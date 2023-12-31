import os
from pathlib import Path

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_file


def does_aizynthfinder_filter_exist() -> bool:
    if not os.path.exists(f"{path_to_data_folder}/aizynthfinder/filter_policy_all.hdf5"):
        return False
    return True

def download_aizynthfinder_filter():
    aizynthfinder_filter = "https://figshare.com/ndownloader/files/25584743"

    directory = f"{path_to_data_folder}/aizynthfinder"
    Path(directory).mkdir(parents=True, exist_ok=True)

    filename = "filter_policy_all.hdf5"
    filepath = f"{directory}/{filename}"
    download_file(aizynthfinder_filter, filepath)