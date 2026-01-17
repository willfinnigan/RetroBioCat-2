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

        self.model = load_model(filename, custom_objects=CUSTOM_OBJECTS, compile=False)

        # Detect model type for Keras 3 compatibility
        self._is_sequential = self.model.__class__.__name__ == 'Sequential'

        # Get dimensions from layer weights (skip layers without weights like InputLayer)
        first_layer_with_weights = None
        last_layer_with_weights = None
        for layer in self.model.layers:
            weights = layer.get_weights()
            if weights:
                if first_layer_with_weights is None:
                    first_layer_with_weights = weights[0]
                last_layer_with_weights = weights[0]

        self._model_dimensions = int(first_layer_with_weights.shape[0]) if first_layer_with_weights is not None else 2048
        self.output_size = int(last_layer_with_weights.shape[1]) if last_layer_with_weights is not None else 1

    def __len__(self):
        return self._model_dimensions

    def predict(self, *args: np.ndarray, **kwargs: np.ndarray):
        # Handle multi-input models (functional API)
        if len(args) > 1 or kwargs:
            # Multi-input model - pass inputs directly
            if kwargs:
                result = self.model.predict(kwargs, verbose=0)
            else:
                result = self.model.predict(list(args), verbose=0)
            return result

        # Single input - may need dimension adjustment for Sequential models saved with Keras 2
        inputs = args[0]
        if self._is_sequential and isinstance(inputs, np.ndarray) and inputs.ndim == 2:
            inputs = np.expand_dims(inputs, axis=1)
        result = self.model.predict(inputs, verbose=0)
        # Squeeze back to 2D if we added a dimension
        if result.ndim == 3 and result.shape[1] == 1:
            result = np.squeeze(result, axis=1)
        return result

