#!python3
##
# @file ToMA.py
# @brief Streamlit GUI entry point for Topological Material Analysis (ToMA).
#
# This is the main entry point for the Streamlit-based graphical user interface.
# Run with: @code{.sh} python -m streamlit run src/ToMA.py @endcode
#
# The GUI provides multi-page navigation for:
# - Single configuration analysis
# - Batch mode processing
# - Circular max flow analysis
# - Cobordism analysis
#
# @authors Yossi Bokor Bleile
# @version 1.3.0
# @date March 2026
# @copyright GPL
#
# Copyright (c) 2023, 2024, 2025 Yossi Bokor Bleile. All rights reserved.

import streamlit as st
# import oineus

st.header("ToMA: Topological Material Analysis")
st.markdown("ToMA is designed to help researchers analyse the structure of materials using tools from geometry and topology. Currently, it uses [Oineus](https://github.com/anigmetov/oineus) and [Diode](https://github.com/mrzv/diode) under the hood to take care of some calculations and processing of the data. As such, it is currently licensed under a GNU GPL license, which you can view from the [License](/License) page in the side bar.")


# welcome_page = st.Page("welcome.py", title="ToMA")
# license_page = st.Page("license.py", title="License")
# single_mode = st.Page("single_mode.py", title="Single Mode")
# multi_mode = st.Page("multi_mode.py", title="Multi Mode")
# batch_mode = st.Page("batch_mode.py", title="Batch Mode")

# pg = st.navigation([welcome_page, single_mode, multi_mode, batch_mode, license_page])#[st.Page("welcome.py", title="ToMA"), st.Page("single_mode.py", title="Single Mode"), st.Page("license.py", title="License")])
# pg.run()