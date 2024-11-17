#!python3
##
# @mainpage Topological Material Analysis (ToMA)
# @authors Yossi Bokor Bleile
# @version 0.4
# @date February 2024
# @copyright GPL
# 
# @section ToMA
# Topological Material Analysis (ToMA) is released under a GPL license. You should have received a [copy](LICENSE.md) of this when you downloaded this repository.

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



st.sidebar.success("Select a mode above.")

st.header("ToMA: Topological Material Analysis")
st.markdown("ToMA is designed to help researchers analyse the structure of materials using tools from geometry and topology. Currently, it uses (Oineus)[https://github.com/anigmetov/oineus] and (Diode) under the hood to take care of some calculations and processing of the data. As such, it is currently licensed under a GNU GPL license, which you can view from the 'License' page in the side bar.")


# welcome_page = st.Page("welcome.py", title="ToMA")
# license_page = st.Page("license.py", title="License")
# single_mode = st.Page("single_mode.py", title="Single Mode")
# multi_mode = st.Page("multi_mode.py", title="Multi Mode")
# batch_mode = st.Page("batch_mode.py", title="Batch Mode")

# pg = st.navigation([welcome_page, single_mode, multi_mode, batch_mode, license_page])#[st.Page("welcome.py", title="ToMA"), st.Page("single_mode.py", title="Single Mode"), st.Page("license.py", title="License")])
# pg.run()