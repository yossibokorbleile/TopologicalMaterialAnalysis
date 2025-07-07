##
# @internal
# @file Batch_Mode.py
# @brief Streamlit page for analysing a single configuration.
# @version 0.1
# @date December 2024
# @author: Yossi Bokor Bleile

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

from toma_functions import *

st.header("Batch Configuration Mode")
comp_tab, plot_tab, vis_tab = st.tabs(["Computation", "Plots", "Visuatlisation"]) #create tabs for the various parts
st.session_state.mode = "multi"
st.markdown("Currently under development")
