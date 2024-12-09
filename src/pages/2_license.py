##
# @internal
# @file 2_License.py
# @brief Streamlit page for displaying the license information.
# @version 0.5
# @date December 2024
# @author Yossi Bokor Bleile

import streamlit as st

import os


def license():
	"""! obtain license information and print it
	"""
	license_text = ""
	l_path = ""
	for p in __file__.split("/")[1:-3]:
		l_path = l_path+"/"+p
	print(l_path)
	l_path = l_path + "/LICENSE.md"
	with open(l_path, 'r') as license:
		license_text = license_text +license.read()
	return license_text
	 
def name_bold():
	"""! return the name in bold
	"""
	return "\033[1m"+"Topological Amorphous Material Analysis"+"\033[0m" 

def copyright():
	"""! return copyright information
	"""	
	return "Copyright (C) Yossi Bokor Bleile\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it under certain conditions.\nTo see the liencese conditions, run `./TAMA.py -l`."
	
def help():
	"""! return instructions on how to get help
	"""	
	return "For help run `./TAMA.py -h`."

st.text(license())