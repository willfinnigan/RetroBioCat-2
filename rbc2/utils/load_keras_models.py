""" This code has been modified from the AIZynthfinder repository at https://github.com/MolecularAI/aizynthfinder """
import numpy as np
import functools
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import logging
import os


def tensorflow_imports():
    try:
        from tensorflow.keras.metrics import top_k_categorical_accuracy
        from tensorflow.keras.models import load_model
        import tensorflow
    except:
        raise Exception("Tensorflow not installed")


    top10_acc = functools.partial(top_k_categorical_accuracy, k=10)
    top10_acc.__name__ = "top10_acc"

    top50_acc = functools.partial(top_k_categorical_accuracy, k=50)
    top50_acc.__name__ = "top50_acc"

    CUSTOM_OBJECTS = {"top10_acc": top10_acc, "top50_acc": top50_acc}

    tf_logger = tensorflow.get_logger()
    tf_logger.setLevel(logging.WARNING)
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
    tensorflow.get_logger().setLevel('WARNING')
    try:
        tensorflow.keras.utils.disable_interactive_logging()
    except:
        pass

    return load_model, CUSTOM_OBJECTS


class LocalKerasModel:
    def __init__(self, filename):
        load_model, CUSTOM_OBJECTS = tensorflow_imports()

        self.model = load_model(filename, custom_objects=CUSTOM_OBJECTS)
        try:
            self._model_dimensions = int(self.model.input.shape[1])
        except AttributeError:
            self._model_dimensions = int(self.model.input[0].shape[1])
        self.output_size = int(self.model.output.shape[1])

    def __len__(self):
        return self._model_dimensions

    def predict(self,  *args: np.ndarray, **_: np.ndarray):
        result = self.model(args)  # this is by far the fastest
        #result = self.model.predict(args)
        #gc.collect()
        return result

