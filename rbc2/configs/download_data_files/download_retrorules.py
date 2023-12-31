import os

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_gdrive_file


def does_retrorules_db_exist():
    if not os.path.exists(f"{path_to_data_folder}/retrorules/retrorules.db"):
        return False
    return True

def download_retrorules_db():
    directory = f"{path_to_data_folder}/retrorules"
    os.makedirs(directory, exist_ok=True)

    file_id = "1CbTQZk-m-VKImMFseXz6kECuDjDXKNew"
    filename = "rr_smi.h5"
    filepath = f"{directory}/{filename}"
    download_gdrive_file(file_id, filepath)

    file_id = "1C7Eb7tiFB1lfDPlqpen8rp6r0pPRf84a"
    filename = "retrorules.db"
    filepath = f"{directory}/{filename}"
    download_gdrive_file(file_id, filepath)

if __name__ == '__main__':
    download_retrorules_db()