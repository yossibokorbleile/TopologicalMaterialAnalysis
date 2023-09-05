# Topological Amorphous Material Analysis

Quick overview: there is a python script Rasmus provided a while ago which does the persistent homology of the ZIF structures. This GUI is meant to first implement this script, and then add certain functionalities. Currently, you need the following Python packages to be installed:
* PySimpleGUI
* ASE
* Matplotlib
* Pandas
* Numpy
* SciPy
* Dionysus
* Diode
* Color
* Oineus
* + any i have forgotten

You will also need to ensure tkinter, cgal, tbb, +others? are installed. 

The GUI currently has two modes: `-s` for looking at a single configureation at a specified time, and `-m` for looking at a single structure at multiple time samples. Functionality for looking at multiple structures at multiple time steps will be added, in particular one any bugs have been ironed out here, and other features added. Or once I come up with a good way of doing it without having to duplicate future work.

## Running the GUI:
Once you have all the packages installed and are in your python environment of choice, navigate to the GUI directory, and from a terminal (on Mac/Linux) or Command prompt(?) (Windows), run 
```python
src/TAMA.py -m s
```
for single sample mode or 

```python
src/TAMA.py -m m
```
for multisample mode or 
```python
src/TAMA.py -m b
```
for batch mode.	
Using the GUI should be pretyy self-explanatory from here, but I can write more things if you would like.

## Documentation
The documentation is built using [doxygen](https://www.doxygen.nl/), which you will need to install if you want to build the documentation locally. Once installed, you just need to run 
```shell
doxygen
```
in the top-level repository directory. Then, the HTML documents will be generated in `Documentation/html`. [doxygen](https://www.doxygen.nl/) also builds LaTeX documentation, which is in `Documentation/latex` which needs to be compiled. Do to so, run 

```
make
``` 
in `Documentation/latex`.


# TODO
* finish visualisation of a loop
* save the outputs in a way which can be relaoded for analysis
* move away from CGAL? (understand alpha complexes and do periodically?)
* include kernel/image/corkernl
* do everything via Oineus?
  
