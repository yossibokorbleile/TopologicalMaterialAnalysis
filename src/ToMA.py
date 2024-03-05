#!python3
##
# @mainpage Topological Material Analysis (ToMA)
# @authors Yossi Bokor Bleile
# @version 0.3
# @date February 2024
# @copyright GPL
# 
# [![DOI](https://zenodo.org/badge/682051112.svg)](https://zenodo.org/doi/10.5281/zenodo.10781424)
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


import argparse
import info
import configparser


parser = argparse.ArgumentParser(prog="Topological Material Analysis", 
									formatter_class = argparse.RawTextHelpFormatter,       
									usage="""This program can be used to analyse atomic figurations to understand porosity of the structures. The main tool used is persistent homology, and its varios related representations. 
									""")
parser.add_argument("-l", "--license", action="store_true", help="""Display the license""")
parser.add_argument("-i", "--interface", choices=["g", "c"], default="g", help="""Select which user interface you would like to use:\n g	use the graphical user interface (Default)\n c	use the command line interface""")
parser.add_argument("-s", "--settings", action="store", help="""Specify the file which contains the settings, only used with the command line interface.""")
parser.add_argument("-n", "--name", action="store", help="""Specify name of the settings to be used, only used with the command line interface.""")
args = parser.parse_args()

print(info.name_bold())
print(info.copyright())
print(info.help())
if args.license:
	print(info.license())
elif args.interface == "g":
	if args.settings != None:
		print("Settings file ignored when using the GUI.")
	from gui import entry_window
	entry_window()
elif args.interface == "c":
	if args.settings == None:
		print("Please specify a file with the settings to use.")
	elif args.name == None:
		print("Please specify a the name of the settings to use.")
	else:	
		mode_config = configparser.ConfigParser()
		mode_config.read()


