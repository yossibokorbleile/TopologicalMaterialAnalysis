# README

[![DOI](https://zenodo.org/badge/682051112.svg)](https://zenodo.org/doi/10.5281/zenodo.10781424)

TopologicalMaterialAnalysis is a Python application to analyse topological structures in materials. It can be used via a graphical user interface or a command line interface. It was developed in conjunction with the [Glass Structure and Mechanics Group, Aalborg Unviersity](https://sites.google.com/view/smedskjaer).

Currently, you need the following Python packages to be installed:
* Streamlit
* ASE
* Pandas
* Numpy
* SciPy
* Dionysus
* Diode
* Oineus
* Plotly Express
* colour
* + any i have forgotten

You will also need to ensure tkinter, CGAL, oneTBB are installed. 

## Documentation
The full documentation is built using [doxygen](https://www.doxygen.nl/), which you will need to install if you want to build the documentation locally. This README contains a quick start guide, and the documentation shouldn't be necessary unless you encounter erros, or wish to contribute to the application.


To build the documentation locally, run 
```shell
doxygen
```
in the top-level repository directory. Then, the HTML documents will be generated in `Documentation/html`. [doxygen](https://www.doxygen.nl/) also builds LaTeX documentation, which is in `Documentation/latex` which needs to be compiled. Do to so, run 

```
make
``` 
in `Documentation/latex`.

### Running the GUI:
Once you have all the packages installed and are in your python environment of choice, and from a terminal (on Mac/Linux) or Command prompt(?) (Windows), run 
```python
streamlit run src/ToMA.py
```

You can then select from 3 different modes: `single`, `multi` and `batch`.
Using the GUI should be pretyy self-explanatory from here.


### Running the CLI: (deprecated)
To run the CLI, use
```python
src/ToMA.py -i c -s $SETTINGSFILE -n $SETTINGSNAME
```

The `SETTINGSFILE` should be an `.ini` file, and each `section` corresponds to a set of variables. Depending on what mode is being used, the contents of the section will differ a little bit. Each `section` needs to contain a variable `MODE` which can have values `single`, `multi`, or `batch`. All three modes also need information about the structure settings to use, which are also stored in an `.ini` file. This structure file can also be used to load structures when using the GUI. For information on formating the structure file, see Structure files section.


### Settings file
To use the CLI of ToMA, it is easiest to specify the settings in a `.ini` file. Depending on which mode of ToMA you want to use, there are different parameters required. Common to them all is `MODE`, which, unsurprisingly, specifies the mode.

## Settings for single mode
In single mode, ToMA analyses a single structure at a single timestep, and so requires the following parameters:
- `STRUCTURE_FILE`: this is the file containing the structure you want to analyse.
- `FILE_FORMAT`: format of the structure file.
- `SAMPLE_TIME`: which time step you want to analyse.
- `CONFIGURATION_FILE`: the file which contains information about configuration you are analysing, see the Configuration files section.
- `CONFIGURATION_NAME`: name of the configuration in the configuration file.
while the following are optional:
- `N_THREADS`: number of threads to use for Oineus.
- `KERNEL`: boolean, set to true if you want to compute kernel persistence.
- `IMAGE`: boolean, set to true if you want to compute image persistence.
- `COKERNEL`: boolean, set to true if you want to compute cokernel persistence.
- `THICKNESS`: set how thick the top and bottom slices should be for kernel/image/cokernel persistence, as a decimal, representing the relative thickness.
- `SAVE_PLOTS`: boolean, set to true if you want to save the plots.

For example:
```
[SINGLE]
MODE = SINGLE
STRUCTURE_FILE = examples/ZIF_example.xyz
SAMPLE_TIME = 0
CONFIGURATION_FILE = examples/structure-types.ini
CONFIGURATION_NAME = ZIF
N_THREADS=32
FILE_FORMAT = xyz
KERNEL = TRUE
IMAGE = TRUE
COKERNEL = TRUE
SAVE_PLOTS = TRUE
THICKNESS = 0.1
N_THREADS = 16
```


## Settings for multi mode
In mutli mode, ToMA analyses a single structure at multiple timesteps, and so requires the following parameters:
- `STRUCTURE_FILE`: this is the file containing the structure you want to analyse.
- `FILE_FORMAT`: format of the structure file.
- `SAMPLE_START`: which time step you begin analysing at.
- `SAMPLE_END`: which time step to end at.
- `SAMPLE_STEPS`: step size.
- `CONFIGURATION_FILE`: the file which contains information about configuration you are analysing, see the Configuration files section.
- `CONFIGURATION_NAME`: name of the configuration in the configuration file.
while the following are optional:
- `N_THREADS`: number of threads to use for Oineus.
- `KERNEL`: boolean, set to true if you want to compute kernel persistence.
- `IMAGE`: boolean, set to true if you want to compute image persistence.
- `COKERNEL`: boolean, set to true if you want to compute cokernel persistence.
- `THICKNESS`: set how thick the top and bottom slices should be for kernel/image/cokernel persistence, as a decimal, representing the relative thickness.
- `SAVE_PLOTS`: boolean, set to true if you want to save the plots.

For example:
```
[MULTI]
MODE = MULTI
STRUCTURE_FILE = examples/ZIF_example.xyz
STRUCTURE_NAME = ZIF
SAMPLE_START = 0
SAMPLE_END = 10
SAMPLE_STEP = 1
CONFIGURATION_FILE =  examples/structure-types.ini
FORMAT = xyz
```

## Settings for batch mode
In batch mode, ToMA analyses every structure it can find in a directory and all subdirectories, at multiple timesteps, and so requires the following parameters:
- `PARENT_DIR`: which directory to commence the search in.
- `FILE_EXT`: what the file extension is, any file with a different extension will be ignored.
- `FILE_FORMAT`: format of the structure file.
- `SAMPLE_START`: which time step you begin analysing at.
- `SAMPLE_END`: which time step to end at.
- `SAMPLE_STEPS`: step size.
- `CONFIGURATION_FILE`: the file which contains information about configuration you are analysing, see the Configuration files section.
- `CONFIGURATION_NAME`: name of the configuration in the configuration file.
while the following are optional:
- `N_THREADS`: number of threads to use for Oineus.
- `KERNEL`: boolean, set to true if you want to compute kernel persistence.
- `IMAGE`: boolean, set to true if you want to compute image persistence.
- `COKERNEL`: boolean, set to true if you want to compute cokernel persistence.
- `THICKNESS`: set how thick the top and bottom slices should be for kernel/image/cokernel persistence, as a decimal, representing the relative thickness.
- `SAVE_PLOTS`: boolean, set to true if you want to save the plots.

For example:
```
[BATCH]
MODE = BATCH
PARENT_DIR = Series_2
FILE_EXT = .transformed
FILE_FORMAT = lammps-dump-text
SAMPLE_START = 100
SAMPLE_END = 200
SAMPLE_STEP = 8
CONFIGURATION_FILE =  structure.ini
CONFIGURATION_NAME = LI2S-P2S5
KERNEL = TRUE
IMAGE = TRUE
COKERNEL = TRUE
SAVE_PLOTS = TRUE
THICKNESS = 0.1
N_THREADS = 16
```

### Configuration files
Configuration files contain information about the type of material you want to analyse. They should be `.ini` files. In particular, it needs to contain information about the atoms you want to consider, the radii to use for each one, and settings for repeating the unit cell. In particular, if you have simulated a material through which certain atoms flow, and want to consider the *backbone* through which these atoms move, you can do not need to modify the file from the simulation: ToMA will only use atoms of the specified types when constructing the $\alpha$-complexes for persistent homology. 

Very important, a structure file can contain information about several different structures, they just need unique names. An example of specifying a structure is:
```
[ZIF]
ATOMS = H, C, N, Zn		
RADII = 0.389, 0.718, 0.635, 1.491
REPEAT_X = 3
REPEAT_Y = 3
REPEAT_Z = 3
```

