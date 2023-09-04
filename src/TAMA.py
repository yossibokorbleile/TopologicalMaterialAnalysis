#!python3
##
# @mainpage Topological Amorphous Material Analysis (AMA)
# @authors Yossi Bokor Bleile
# @version 0.1
# @date April 2023
# @copyright if CGAL still used, GPL
#
# @section description_main Description
# Topological Amorphous Material Analysis is relies on Dionysus, Oineus and (currently) Diode to handle a lot of the "under the hood" functionalities.
# generating source code doca#!/usr/bin/python3
##
# @mainpage Amorphous Material Analysis (AMA)
# @authors Yossi Bokor Bleile
# @version 0.1
# @date April 2023
# @copyright if CGAL still used, GPL
#
# @section description_main Description
# Topological Amorphous Material Analysis is relies on Dionysus, Oineus and (currently) Diode to handle a lot of the "under the hood" functionalities.
# generating source code documentation with Doxygen.
#
# @section arg_parser Argument Parser for Amorphous Material Analysis
# There are two ways you can use Amorphous Material Analysis:
# -# Graphical User Interface (default)
# -# Command Line Interface
#    
# both of which have several different modes:a
# -# single (default)
# -# multi
# -# batch
#
# which can be selected when starting AMA using flags and options:
# - "-i"    interface flag
#   - "g"   graphical user interface option
#   - "c"   command line interface option
# - "-m"    mode flag
#   - "s"   single mode option
#   - "m"   multi mode option
#   - "b"   batch mode option
#
# To see which modes are availble in the two interfaces, see gui.py and cli.py respectively,
#
# Copyright (c) 2023 Yossi Bokor Bleile.  All rights reserved.    


import argparse
import info


parser = argparse.ArgumentParser(prog="Topological Amorphous Material Analysis", 
                                    formatter_class = argparse.RawTextHelpFormatter,       
                                    usage="""This program can be used to analyse atomic figurations to understand porosity of the structures. The main tool used is persistent homology, and its varios related representations. 
                                    """)
parser.add_argument("-l", "--license", action="store_true", help="""Display the license""")
parser.add_argument("-i", "--interface", choices=["g", "c"], default="g", help="""Select which user interface you would like to use:\n g	use the graphical user interface (Default)\n c	use the command line interface""")
parser.add_argument("-m", "--mode", choices=["s", "m", "b"], default="s", help="""Select which mode you would like to launch\n s	single analysis mode: use this mode to analyse a single initial atom confiugration, at a speicifc time (Default)\n m	multi analysis mode: use this mode to analyse a single initial atom configuration, at several time steps\n b	batch analysis mode: use this mode to analyse several initial atom configations""")
parser.add_argument("-d", "--directory", help="""Specify the directory which contains the initial configuration files for batch mode. ONLY BATCH MODE""")
parser.add_argument("-e", "--extension", action="store", help="Specify the extension the initial configuration files have. ONLY BATCH MODE")
parser.add_argument("-f", "--file", action="store", help="""Specify the initial configuration file. ONLY IN MULTI AND SINGLE MODE""")
parser.add_argument("--format", action="store", default="xyz", help="""Specify the format of the configuration file.""")
parser.add_argument("-s", "--settings", action="store", help="""Specify the file which contains the settings.""")
parser.add_argument("-t", "--type", help="""Select the structure type you are analysing.""")
parser.add_argument("-r", "--range", action="store", help="""Specify the sample range settings to use.""")
parser.add_argument("-k", "--kernel", action="store", help="""Specify kernel/image/cokernel settings.""")
parser.add_argument("--sample-time", type=int, action="store", help="""Time at which to sample the configuratio.""")
parser.add_argument("--save-all", action="store_true",dest="save_all", help="""Save the images for each sample time. ONLY BATCH AND MULTI MODE""")
args = parser.parse_args()

info.intro()

if args.license:
	info.license()
 
if args.mode == "s" and args.interface == "g":
    from gui import single_mode as mode
    mode()
elif args.mode == "m" and args.interface == "g":
    from gui import multi_mode as mode
    mode()
elif args.mode == "b" and args.interface == "g":
    from gui import batch_mode as mode
    mode()
elif args.mode == "s" and args.kernel == None and args.interface == "c":
    from cli import single_mode as mode
    mode(args.file, args.format, args.settings, args.type,  args.sample_time)
elif args.mode == "s" and args.kernel != None and args.interface == "c":
	from cli import single_mode_kernel as mode
	mode(args.file, args.settings, args.format, args.type, args.sample_time, args.kernel)
elif args.mode == "m" and args.kernel == None and args.interface == "c":
    from cli import multi_mode as mode
    mode(args.file, args.settings, args.format, args.type, args.range, save_all=args.save_all)
elif args.mode == "m" and args.kernel != None and args.interface == "c":
    from cli import multi_mode_kernel as mode
    mode(args.file, args.settings, args.format, args.type, args.range, args.kernel, save_all=args.save_all)
elif args.mode == "b" and args.kernel == None and args.interface == "c":
    from cli import batch_mode as mode 
    mode(args.directory, args.extension, args.format, args.settings, args.type, args.range, save_all=args.save_all)
elif args.mode == "b" and args.kernel != None and args.interface == "c":
	from cli import batch_mode_kernel as mode
	mode(args.directory, args.extension, args.format, args.settings, args.type, args.range, args.kernel, save_all=args.save_all)

