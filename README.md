# README

[![DOI](https://zenodo.org/badge/682051112.svg)](https://zenodo.org/doi/10.5281/zenodo.10781424)

TopologicalMaterialAnalysis is a Python application to analyse topological structures in materials. It can be used via a graphical user interface or a command line interface. It was developed in conjunction with the [Glass Structure and Mechanics Group, Aalborg Unviersity](https://sites.google.com/view/smedskjaer).

Currently, you need the following Python packages to be installed:
* PySimpleGUI
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
Once you have all the packages installed and are in your python environment of choice, navigate to the GUI directory, and from a terminal (on Mac/Linux) or Command prompt(?) (Windows), run 
```python
src/ToMA.py
```

You can then select from 3 different modes: `single`, `multi` and `batch`.
Using the GUI should be pretyy self-explanatory from here.


### Running the CLI:
To run the CLI, use
```python
src/ToMA.py -i c -s $SETTINGSFILE -n $SETTINGSNAME
```

The `SETTINGSFILE` should be an `.ini` file, and each `section` corresponds to a set of variables. Depending on what mode is being used, the contents of the section will differ a little bit. Each `section` needs to contain a variable `MODE` which can have values `single`, `multi`, or `batch`. All three modes also need information about the structure settings to use, which are also stored in an `.ini` file. This structure file can also be used to load structures when using the GUI. For information on formating the structure file, see [Structure files](#Structure-files)

#### Settings for single mode

### Structure files


