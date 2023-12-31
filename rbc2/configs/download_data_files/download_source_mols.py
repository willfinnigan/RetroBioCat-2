import os

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_file


def does_source_mols_db_exist():
    if not os.path.exists(f"{path_to_data_folder}/buyability/source_mols.db"):
        return False
    return True

def download_source_mols_db():
    source_mols_db = "https://figshare.com/ndownloader/files/43860240"

    directory = f"{path_to_data_folder}/buyability"
    os.makedirs(directory, exist_ok=True)

    filename = "source_mols.db"
    filepath = f"{directory}/{filename}"
    download_file(source_mols_db, filepath)

if __name__ == '__main__':
    download_source_mols_db()
