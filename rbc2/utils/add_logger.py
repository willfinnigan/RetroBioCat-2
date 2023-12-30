import logging
import os

handler = logging.StreamHandler()
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

def add_logger(name, level='DEBUG'):
    logger = logging.getLogger(name)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False
    return logger


def set_tensorflow_logging(level='2'):
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = level