import gzip
import os
import shutil
from pathlib import Path

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_file


def does_askcos_exist() -> bool:
    if not os.path.exists(f"{path_to_data_folder}/askcos/reaxys/retro.templates.json.gz"):
        return False
    if not os.path.exists(f"{path_to_data_folder}/askcos/reaxys/variables/variables.data-00000-of-00001"):
        return False
    return True

def download_askcos_model():
    askcos_model = "https://github.com/ASKCOS/askcos-data/blob/main/models/template_prioritization/reaxys/1/variables/variables.data-00000-of-00001?raw=true"
    askcos_templates = "https://github.com/ASKCOS/askcos-data/blob/main/templates/retro.templates.json.gz?raw=true"

    # if askcos folder doesn't exist, create it with Pathlib
    directory = f"{path_to_data_folder}/askcos/reaxys"

    Path(directory).mkdir(parents=True, exist_ok=True)

    filename = "variables.data-00000-of-00001"
    filepath = f"{directory}/variables/{filename}"
    download_file(askcos_model, filepath)

    filename = "retro.templates.json.gz"
    filepath = f"{directory}/{filename}"
    download_file(askcos_templates, filepath)


if __name__ == '__main__':
    download_askcos_model()