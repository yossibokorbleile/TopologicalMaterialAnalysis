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
		mode_config.read(args.settings)
		mode = mode_config.get(args.name, "MODE")
	if mode == "SINGLE":
		from cli import single_mode
		structure_file = mode_config.get(args.name,"STRUCTURE_FILE")
		file_format = mode_config.get(args.name,"FILE_FORMAT")
		configuration_file = mode_config.get(args.name,"CONFIGURATION_FILE")
		configuration = mode_config.get(args.name,"CONFIGURATION_NAME")
		sample_time = mode_config.get(args.name,"SAMPLE_TIME")
		try:
			n_threads=mode_config.get(args.name,"N_THREADS")
		except:
			n_threads = 4
		try:
			save_plots = mode_config.get(args.name,"SAVE_PLOTS")
		except:
			save_plots = False
		try:
			kernel = mode_config.get(args.name,"KERNEL")
		except:
			kernel = False
		try:
			image = mode_config.get(args.name,"IMAGE")
		except:
			image = False
		try:
			cokernel = mode_config.get(args.name,"COKERNEL")
		except:
			cokernel = False
		try:
			upper_threshold = mode_config.get(args.name, "UPPER_THRESHOLD")
		except:
			upper_threshold = "Auto"
		try:
			lower_threshold = mode_conifg.get(args.name, "LOWER_THRESHOLD")
		except:
			lower_threshold = "Auto"
		single_mode(structure_file=structure_file, file_format=file_format,configuration_file=configuration_file, configuration=configuration, sample_time=sample_time,  n_threads=n_threads, save_plots=save_plots, kernel=kernel, image=image, cokernel=cokernel,upper_threshold=upper_threshold, lower_threshold=lower_threshold)
	elif mode == "MULTI":
		from cli import multi_mode
	elif mode == "BATCH":
		from cli import batch_mode


	# 	config.read(configuration_file) #load the file containing the structures
	# atoms = [str(a).strip() for a in config.get(configuration, "ATOMS").split(",")]



