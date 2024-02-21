#!python3
##
# @mainpage Topological Amorphous Material Analysis (TAMA)
# @authors Yossi Bokor Bleile
# @version 0.3
# @date February 2024
# @copyright GPL
# 
# @section License
# Topological Amorphous Material Analysis (TAMA) is released under a GPL license. You should have received a [copy](LICENSE.md) of this when you downloaded this repository.
# @section description_main Description
# Topological Amorphous Material Analysis is relies on Dionysus, Oineus and (currently) Diode to handle a lot of the "under the hood" functionalities.
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
# In the GUI, these are selected at startup, and for the CLI, the mode to be used is set in a settings file with an INI format. To use `TAMA` in the CLI mode, run 
#`./TAMA.py -i c -s $SETTINGS FILE`.
#
# Information about the packages required and instructions for the settings files are in the [README](README.md).
#
# Copyright (c) 2023, 2024 Yossi Bokor Bleile.  All rights reserved.    


import argparse
import info
import configparser


parser = argparse.ArgumentParser(prog="Topological Amorphous Material Analysis", 
									formatter_class = argparse.RawTextHelpFormatter,       
									usage="""This program can be used to analyse atomic figurations to understand porosity of the structures. The main tool used is persistent homology, and its varios related representations. 
									""")
parser.add_argument("-l", "--license", action="store_true", help="""Display the license""")
parser.add_argument("-i", "--interface", choices=["g", "c"], default="g", help="""Select which user interface you would like to use:\n g	use the graphical user interface (Default)\n c	use the command line interface""")
parser.add_argument("-s", "--settings", action="store", help="""Specify the file which contains the settings, only used with the command line interface.""")
parser.add_argument("-n", "--name", action="store", help="""Specify name of the settings to be used, only used with the command line interface.""")
args = parser.parse_args()
print(args.settings)
info.intro()

if args.license:
	info.license()
 
if args.interface == "g":
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


