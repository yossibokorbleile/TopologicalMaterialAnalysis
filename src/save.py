##
# @internal
# @file save.py
# @brief Functions to save various outputs. 

import plotly.express as px
import os

def save_all_single(dir : str, config_name : str, sample : int, births : list, deaths : list, kicr, values_main, ):