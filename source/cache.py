#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/cache.py

# List general Python imports
import sys
import os
import pickle


def save_object(obj, name, path='/dev/shm/addtag'):
    if isinstance(name, str):
        if (len(name) > 0):
            os.makedirs(path, exist_ok=True)
            with open(os.path.join(path, name), 'wb') as flo:
                pickle.dump(obj, flo, pickle.HIGHEST_PROTOCOL)
        else:
            raise Exception("object name is not a string")

def load_object(name, path='/dev/shm/addtag'):
    with open(os.path.join(path, name), 'rb') as flo:
        data = pickle.load(flo)
    return data

