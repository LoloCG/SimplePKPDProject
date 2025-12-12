C++ CLI pharmakokinetics One compartment (IV/EV) model.
Used along simple python GUI with matplotlib. 

### Why
Reason for creation is mostly learning purposes:
- Learning of C++ with OOP, CMake, and related.
- Practice of basic PK compartmental modeling.
- Understand mathematical implications of flip-flop kinetics, tlag, etc.


### Usage
Compiled c++ .exe accepts basic pk arguments, giving amounts in csv format (CLI printed or file format).
Python GUI script can run .exe directly, described in EXE_PATH global var. PK arguments can be described in ARGS global variable.

EV administration is selected when providing `ka` or with `--ev` argument (defaults to 0.1).
More usage info can be obtained with `--help` or `-h` flags.
