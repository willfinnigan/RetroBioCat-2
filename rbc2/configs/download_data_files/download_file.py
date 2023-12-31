import requests
from tqdm import tqdm
import gdown

def download_file(url, filename):
    """
    Downloads a file from a given URL and saves it with the given filename, with a progress bar.

    Args:
    url: The URL of the file to download.
    filename: The name to save the file as.

    Returns:
    None
    """
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        file_size = int(response.headers.get('Content-Length', 0))
        progress = tqdm(response.iter_content(1024), f'Downloading {filename}', total=file_size, unit='B', unit_scale=True, unit_divisor=1024)
        with open(filename, 'wb') as file:
            for data in progress.iterable:
                file.write(data)
                progress.update(len(data))
    else:
        print("Error: Received response code " + str(response.status_code) + " from server.")




def download_gdrive_file(file_id, output):
    """
    Downloads a file from Google Drive.

    Args:
    file_id: The ID of the file on Google Drive.
    output: The name of the file to save the download as.

    Returns:
    None
    """
    url = f'https://drive.google.com/uc?id={file_id}'
    gdown.download(url, output, quiet=False)