##
# @internal
# @file info.py
# @authors Yossi Bokor Bleile
# @version 0.1
# @date April 2023
# @copyright BSD
#


import os
 
# get current directory
path = os.getcwd()
print("Current Directory", path)
 
# prints parent directory
print(os.path.abspath(os.path.join(path, os.pardir)))

def license():
    """! obtain license information and print it
    """
    with open(os.path.abspath(os.path.join(path, os.pardir))+"/LICENSE", 'r') as license:
    	print(license.read())
     
def intro():
    """! print brief introduction about the program
    """
    print("Topological Material Analysis Copyright (C) Yossi Bokor Bleile\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it under certain conditions.\nTo see the liencese conditions, run `ToMA.py -l`.\nFor help run `ToMA.py -h`.")
    return "Topological Amorphous Material Analysis Copyright (C) Yossi Bokor Bleile, Aalborg University\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it under certain conditions.\nTo see the liencese conditions, run `ToMA.py -l`.\nFor help run `ToMA.py -h`."
	