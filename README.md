Path Integral Molecular Dynamics (PIMD) Simulation Package
==========================================================
Welcome to the Path Integral Molecular Dynamics (PIMD) Simulation Package GitHub repository! This software package offers a comprehensive suite of tools for performing various path integral molecular dynamics simulations, including PIMD (both preconditioned and non-preconditioned), Ring Polymer Molecular Dynamics (RPMD), Thermostatted RPMD, and classical MD simulations. Additionally, it supports umbrella sampling simulations for enhanced sampling of complex systems and contains Fortran routines for use with Plumed software [Tribello *et al.*](https://doi.org/10.1016/j.cpc.2013.09.018)

Features:
=========
  1. **Quantum thermodynamics properties:** Includes preconditioned and non-preconditioned PIMD methods (see [Korol *et al.*](https://doi.org/10.1063/1.5134810), [Ceriotti *et al.*](http://dx.doi.org/10.1063/1.3489925), [Shiga](10.1016/B978-0-12-409547-2.11614-2) ) for accurate simulation of static (thermodynamics) quantum properties in molecular systems. 
  2. **Approximate quantum dynamical properties:** Implements Ring Polymer Molecular Dynamics (RPMD) [Craig & Manolopoulos](https://doi.org/10.1063/1.1777575)  and Thermostatted Ring Polymer Molecular Dynamics (T-RPMD) [Rossi & Manolopoulos](https://doi.org/10.1063/1.4883861) for studying quantum dynamical effects in molecular systems.
  3. **Free Energy Calculations:** Facilitates enhanced sampling simulations using umbrella sampling techniques. Currently supports bonds, angles, dihedrals, and their differences as collective variables. Additional collective variables can be easily added for specialized simulations. It also includes Fortran routines for interfacing with Plumed [Tribello *et al.*](https://doi.org/10.1016/j.cpc.2013.09.018) software for free energy calculation.
  4. **Parallelization:** Utilizes OpenMP for efficient parallelization, enhancing performance on multi-core systems.
  5. **Custom PES:** Allows custom potential energy surfaces (PESs), enabling simulations tailored to specific systems or research questions.
  6. **Available Thermostats:** Currently supports Massive Nose-Hoover chain thermostat (MNHC) [Tuckerman *et al.*](http://dx.doi.org/10.1080/00268979600100761)  and path integral Langevin equation (PILE) thermostat [Ceriotti et al.](http://dx.doi.org/10.1063/1.3489925) for PIMD simulations. The Nose-Hoover chain and the Langevin equation thermostats are available for classical simulations.
  7. **Ab-initio PIMD and Classical MD:** Supports ab-initio path integral and classical MD simulations. Currently, it is only implemented for [Gaussian](https://gaussian.com/) software.

Prerequisites:
==============
It requires Intel MKL libraries. 

Getting Started
===============
To get started with using the PIMD Simulation Package, follow these steps:

  Clone the Repository:
  ---------------------
    git clone https://github.com/arandharamrinal/PIMD.git

  Compile the Code 
  -----------------
    #(make necessary changes in the Makefile).
    cd PIMD  
    make
    #This will generate an executable in the bin directory.
  Run a Simulation:
  -----------------
    #Copy pimd.exe to the directory you want to run the simulation.
    ./pimd.exe input_parameter_file

Explore Documentation:
======================
Detailed documentation describing input file format, available options, and usage examples will be available soon. If you have any questions, please feel free to contact me directly.

Contributing:
============
We welcome contributions from the community to improve and extend the capabilities of the PIMD Simulator. If you have suggestions or bug reports or would like to contribute code, please open an issue or submit a pull request.

Contact
=======
For any questions, assistance, or collaboration opportunities, please don't hesitate to contact the repository owner:

Name: [Mrinal Arandhara]
Email: [arandharamrinal@gmail.com]

License
This software is distributed under the [MIT license]. See the LICENSE file for details.
