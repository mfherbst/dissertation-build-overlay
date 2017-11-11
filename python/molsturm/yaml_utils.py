import numpy as np
import yaml


def construct_ndarray(loader, node):
    """YAML constructor for numpy arrays"""
    value = loader.construct_sequence(node, deep=True)
    return np.array(value)


def install_constructors():
    """Install all YAML constructors defined in this module"""
    yaml.constructor.SafeConstructor.add_constructor("!ndarray", construct_ndarray)
