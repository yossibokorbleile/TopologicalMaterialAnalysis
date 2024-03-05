# README
[![DOI](https://zenodo.org/badge/682051112.svg)](https://zenodo.org/doi/10.5281/zenodo.10781424)
Quick overview: there is a python script Rasmus provided a while ago which does the persistent homology of the ZIF structures. This GUI is meant to first implement this script, and then add certain functionalities. Currently, you need the following Python packages to be installed:
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

You will also need to ensure tkinter, cgal, tbb, +others? are installed. 

The GUI currently has two modes: `-s` for looking at a single configureation at a specified time, and `-m` for looking at a single structure at multiple time samples. Functionality for looking at multiple structures at multiple time steps will be added, in particular one any bugs have been ironed out here, and other features added. Or once I come up with a good way of doing it without having to duplicate future work.

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


