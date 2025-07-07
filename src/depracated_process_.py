##
# @internal
# @file process.py
# @brief functions to process the data using oineus and diode.
# @version 0.5
# @date December 2024

import streamlit as st
# import oineus
import numpy as np
import pandas as pd
import os

import configparser
from ase import io, Atoms
import diode
import math
from colour import Color
from scipy.interpolate import interpn
from functools import cmp_to_key
from scipy.interpolate import interpn
import plotly.express as px
import plotly.graph_objects as go

from toma_io import *


