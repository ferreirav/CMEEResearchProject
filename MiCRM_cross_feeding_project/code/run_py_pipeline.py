#!/usr/bin/env python3

"""
Doc Strings
"""

__appname__ = 'run_py_pipeline.py'
__author__ = 'Vitor Ferreira (f.ferreira22@imperial.ac.uk)'
__version__ = '0.0.1'


# IMPORTS
import os
import time
import sys
import subprocess
import pandas as pd
from multiprocessing import Pool
from emp_data_wrang import get_data, get_genome
import carveme
from smetana.interface import run_detailed