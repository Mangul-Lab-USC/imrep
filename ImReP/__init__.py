
import os

import ImReP.utils

import ImReP.cast
from ImReP.cast import *

import ImReP.settings
from ImReP.settings import *

import ImReP.imrep
from .imrep import *

import ImReP.info

_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_db(path):
    return os.path.join(_ROOT, 'db', path)
