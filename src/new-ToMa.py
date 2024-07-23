#!python3
##
# @mainpage Topological Material Analysis (ToMA)
# @authors Yossi Bokor Bleile
# @version 0.3
# @date February 2024
# @copyright GPL
# 
# @section License
# Topological Material Analysis (ToMA) is released under a GPL license. You should have received a [copy](LICENSE.md) of this when you downloaded this repository.
# @section description_main Description
# ToMA is dependent on Oineus and Diode to handle a lot of the "under the hood" functionalities.
#
# There are two ways you can use Amorphous Material Analysis:
# -# Graphical User Interface (default)
# -# Command Line Interface
#    
# both of which have several different modes:
# -# single
# -# multi
# -# batch
#
# In the GUI, these are selected at startup, and for the CLI, the mode to be used is set in a settings file with an INI format. To use `ToMA` in the CLI mode, run 
#`./ToMA.py -i c -s $SETTINGS FILE`.
#
# Information about the packages required and instructions for the settings files are in the [README](README.md).
#
# Copyright (c) 2023, 2024 Yossi Bokor Bleile.  All rights reserved.    


import streamlit as st

welcome_page = st.Page("welcome.py", title="ToMA")
license_page = st.Page("license.py", title="ToMA License")
single_mode = st.Page("single_mode.py", title="ToMA single mode")

pg = st.navigation([welcome_page, single_mode, license_page])
pg.run()