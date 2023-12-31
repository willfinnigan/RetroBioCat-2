import os
import shutil
from pathlib import Path

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_file


def does_ringbreaker_exist() -> bool:
    if not os.path.exists(f"{path_to_data_folder}/ringbreaker/models"):
        return False
    elif not os.path.exists(f"{path_to_data_folder}/ringbreaker/data"):
        return False
    return True


def download_ringbreaker_model():
    ringbreaker_model = "https://figshare.com/ndownloader/files/22789826"
    ringbreaker_data = "https://figshare.com/ndownloader/files/22789832"

    # if ringbreaker folder doesn't exist, create it with Pathlib
    directory = f"{path_to_data_folder}/ringbreaker"

    Path(directory).mkdir(parents=True, exist_ok=True)

    filename = "ringbreaker_models.zip"
    filepath = f"{directory}/{filename}"
    download_file(ringbreaker_model, filepath)
    shutil.unpack_archive(filepath, directory)
    os.remove(filepath)

    filename = "ringbreaker_data.zip"
    filepath = f"{directory}/{filename}"
    download_file(ringbreaker_data, filepath)
    shutil.unpack_archive(filepath, directory)
    os.remove(filepath)


if __name__ == '__main__':
    download_ringbreaker_model()
