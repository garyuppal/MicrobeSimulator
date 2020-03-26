# MicrobeSimulator
General hybrid agent-based/continuous-field simulator for microbes interacting via secreted compounds. With finite element (using deal.II) and finite difference implementations for continuous fields.

## Installation and usage
The program requires the prior installation of the open source deal.ii library available [here](https://www.dealii.org/). A makefile can then be generated using CMake. Use the command `cmake -DDEAL_II_DIR=path/to/dealii/ .` to generate the makefile. Program is set up with a parameter text file. Information on parameters and configurations is given in the configurations directory.

### To Do
- [ ] extend filter and mixer meshes to 3D
- [ ] add vortex and cylinder geometries

