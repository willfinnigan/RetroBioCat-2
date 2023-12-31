import os

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_gdrive_file


def does_enzymemap_exist():
    if not os.path.exists(f"{path_to_data_folder}/enzymemap/brenda/default_retro-weights.hdf5"):
        return False
    return True

def download_enzymemap():
    directory = f"{path_to_data_folder}/enzymemap/brenda"
    os.makedirs(directory, exist_ok=True)

    file_id = "1CNS8AZ49z3V24zXgIJoqh2EWWGtxNPXM"
    filename = "enzymemap_brenda_metadata.csv"
    filepath = f"{directory}/{filename}"
    download_gdrive_file(file_id, filepath)

    file_id = "1CSE1ZszPvhhrRUXJU4YZC9WJSeqF9pXt"
    filename = "enzymemap_brenda_metadata.csv"
    filepath = f"{directory}/{filename}"
    download_gdrive_file(file_id, filepath)

    file_id = "1CAGlHEp4Cnbl_FgepIdXYRCMHd562WDv"
    filename = "default_retro-weights.hdf5"
    filepath = f"{directory}/{filename}"
    download_gdrive_file(file_id, filepath)



if __name__ == '__main__':
    download_enzymemap()