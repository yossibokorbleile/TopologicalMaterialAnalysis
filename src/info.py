##
# @internal
# @file info.py
# @authors Yossi Bokor Bleile
# @version 0.1
# @date April 2023
# @copyright BSD
#

import os

def license():
	"""! obtain license information and print it
	"""
	license_text = ""
	l_path = ""
	for p in __file__.split("/")[1:-2]:
		l_path = l_path+"/"+p
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
