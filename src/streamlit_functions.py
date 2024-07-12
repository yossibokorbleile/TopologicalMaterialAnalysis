import streamlit as st
import process

def load_configuration_settings():
	st.session_state.atoms, st.session_state.radii, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z = process.read_configuration(st.session_state["config_file"], st.session_state["config_name"])

def load_computation_settings():
	st.session_state.n_threads, st.session_state.save_plots, st.session_state.kernel, image, st.session_state.cokernel, st.session_state.thickness = process.load_computation_settings(st.session_state["comp_file"], st.session_state["comp_name"])

