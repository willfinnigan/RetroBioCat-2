import os
import logging
import zipfile
import tempfile
import glob as _glob
from pathlib import Path

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_file import download_file, download_gdrive_file

logger = logging.getLogger(__name__)

# Graph2SMILES (coleygroup) - Google Drive + GitHub
GRAPH2SMILES_MODEL_GDRIVE_ID = "11IJzTvyEvLclVSDN6gz2wy1Pg2hSOqa0"
GRAPH2SMILES_VOCAB_URL = "https://raw.githubusercontent.com/coleygroup/Graph2SMILES/main/preprocessed/default_vocab_smiles.txt"

# Augmented Transformer (ASKCOS) - Dropbox .mar archive
AUGMENTED_TRANSFORMER_MAR_URL = "https://www.dropbox.com/scl/fi/m36a9k2wpqxv37agujx0j/USPTO_FULL.mar?rlkey=6k8ulce7dcjkw8nn9y07mr7ag&dl=1"


# --- EGraph2RSMILES ---
# No pre-trained weights are publicly available for this model.

def _get_egraph2rsmiles_dir():
    return os.path.join(path_to_data_folder, 'template_free', 'egraph2rsmiles')


def does_egraph2rsmiles_exist() -> bool:
    d = _get_egraph2rsmiles_dir()
    return (os.path.exists(os.path.join(d, 'model.pt')) and
            os.path.exists(os.path.join(d, 'vocab.txt')))


def download_egraph2rsmiles():
    logger.warning("EGraph2RSMILES has no publicly available pre-trained weights. "
                   "Please manually place model.pt and vocab.txt in %s",
                   _get_egraph2rsmiles_dir())


def ensure_egraph2rsmiles_exists():
    if not does_egraph2rsmiles_exist():
        download_egraph2rsmiles()
    if not does_egraph2rsmiles_exist():
        raise FileNotFoundError(
            f"EGraph2RSMILES model files not found in {_get_egraph2rsmiles_dir()}. "
            f"No public pre-trained weights are available for this model.")


# --- Graph2SMILES ---

def _get_graph2smiles_dir():
    return os.path.join(path_to_data_folder, 'template_free', 'graph2smiles')


def does_graph2smiles_exist() -> bool:
    d = _get_graph2smiles_dir()
    return (os.path.exists(os.path.join(d, 'model.pt')) and
            os.path.exists(os.path.join(d, 'vocab.txt')))


def download_graph2smiles():
    """Download Graph2SMILES model weights from Google Drive and vocab from GitHub."""
    d = _get_graph2smiles_dir()
    Path(d).mkdir(parents=True, exist_ok=True)

    model_path = os.path.join(d, 'model.pt')
    vocab_path = os.path.join(d, 'vocab.txt')

    if not os.path.exists(model_path):
        logger.info("Downloading Graph2SMILES checkpoint from Google Drive...")
        download_gdrive_file(GRAPH2SMILES_MODEL_GDRIVE_ID, model_path)

    if not os.path.exists(vocab_path):
        logger.info("Downloading Graph2SMILES vocabulary...")
        download_file(GRAPH2SMILES_VOCAB_URL, vocab_path)


def ensure_graph2smiles_exists():
    if not does_graph2smiles_exist():
        logger.info("Graph2SMILES model files not found, downloading...")
        download_graph2smiles()
    if not does_graph2smiles_exist():
        raise FileNotFoundError(
            f"Graph2SMILES model files not found in {_get_graph2smiles_dir()}. "
            f"Download failed - please check your internet connection.")


# --- Augmented Transformer ---

def _get_augmented_transformer_dir():
    return os.path.join(path_to_data_folder, 'template_free', 'augmented_transformer')


def does_augmented_transformer_exist() -> bool:
    d = _get_augmented_transformer_dir()
    checkpoints = _glob.glob(os.path.join(d, 'model_step_*.pt'))
    return len(checkpoints) > 0


def download_augmented_transformer():
    """Download Augmented Transformer .mar archive from Dropbox and extract checkpoint."""
    d = _get_augmented_transformer_dir()
    Path(d).mkdir(parents=True, exist_ok=True)

    logger.info("Downloading Augmented Transformer model from Dropbox...")

    with tempfile.TemporaryDirectory() as tmpdir:
        mar_path = os.path.join(tmpdir, 'USPTO_FULL.mar')
        download_file(AUGMENTED_TRANSFORMER_MAR_URL, mar_path)

        if not os.path.exists(mar_path) or os.path.getsize(mar_path) == 0:
            raise RuntimeError("Failed to download Augmented Transformer .mar file")

        # .mar files are ZIP archives - extract model checkpoint
        logger.info("Extracting model checkpoint from .mar archive...")
        with zipfile.ZipFile(mar_path, 'r') as zf:
            for name in zf.namelist():
                if name.startswith('model_step_') and name.endswith('.pt'):
                    zf.extract(name, d)
                    logger.info("Extracted %s", name)
                    return

        # If no model_step_*.pt found, extract all .pt files
        with zipfile.ZipFile(mar_path, 'r') as zf:
            pt_files = [n for n in zf.namelist() if n.endswith('.pt')]
            if pt_files:
                for name in pt_files:
                    zf.extract(name, d)
                    logger.info("Extracted %s", name)
            else:
                raise RuntimeError(
                    "No .pt checkpoint found in .mar archive. Contents: "
                    + str(zf.namelist()))


def ensure_augmented_transformer_exists():
    if not does_augmented_transformer_exist():
        logger.info("Augmented Transformer model not found, downloading...")
        download_augmented_transformer()
    if not does_augmented_transformer_exist():
        raise FileNotFoundError(
            f"Augmented Transformer model not found in "
            f"{_get_augmented_transformer_dir()}. "
            f"Download failed - please check your internet connection.")


if __name__ == '__main__':
    print(f"EGraph2RSMILES exists: {does_egraph2rsmiles_exist()}")
    print(f"Graph2SMILES exists: {does_graph2smiles_exist()}")
    print(f"Augmented Transformer exists: {does_augmented_transformer_exist()}")
