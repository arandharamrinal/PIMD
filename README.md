# PIMD
Path Integral Molecular Dynamics (PIMD) Simulation Package
==========================================================
Welcome to the Path Integral Molecular Dynamics (PIMD) Simulation Package GitHub repository! This software package offers a comprehensive suite of tools for performing various path integral molecular dynamics simulations, including PIMD (both preconditioned and non-preconditioned), Ring Polymer Molecular Dynamics (RPMD), Thermostatted RPMD, and classical MD simulations. Additionally, it supports umbrella sampling simulations for enhanced sampling of complex systems.

Features:
=========
  *Quantum thermodynamics properties:* Includes preconditioned and non-preconditioned PIMD methods for accurate simulation of static (thermodynamics) quantum properties in molecular systems. 
  *Approximate quantum dynamical properties:* Implements Ring Polymer Molecular Dynamics (RPMD) and Thermostatted RPMD for studying quantum dynamical effects in molecular systems.
  *Free Energy Calculations:* Facilitates enhanced sampling simulations using umbrella sampling techniques. Currently supports bonds, angles, dihedrals, and their differences as collective variables. Additional collective variables can be easily added for specialized simulations. It also includes Fortran routines for interfacing with Plumed software for free energy calculation.
  *Parallelization:* Utilizes OpenMP for efficient parallelization, enhancing performance on multi-core systems.
  *Custom PES:* Allows custom potential energy surfaces (PESs), enabling simulations tailored to specific systems or research questions.
  *Available Thermostats:* Currently supports Massive Nose-Hoover chain thermostat (MNHC)(Tuckerman et al.) and path integral Langevin equation (PILE) thermostat (Ceriotti et al.) for PIMD simulations. The Nose-Hoover chain thermostat and the Langevin equation thermostats are available for classical simulations.


Getting Started
===============
To get started with using the PIMD Simulator, follow these steps:

  Clone the Repository:
  ---------------------
    bash
    git clone https://github.com/arandharamrinal/PIMD.git

  Compile the Code:
  -----------------
    bash
    cd PIMD
    make
  Run a Simulation:
  -----------------
    bash
    ./pimd_simulation input_file.inp

Explore Documentation:
======================
Detailed documentation describing input file format, available options, and usage examples will be available soon. If you have any questions, please feel free to contact me directly.
Contributing
We welcome contributions from the community to improve and extend the capabilities of the PIMD Simulator. If you have suggestions or bug reports or would like to contribute code, please feel free to open an issue or submit a pull request.

Contact
=======
For any questions, assistance, or collaboration opportunities, please don't hesitate to contact the repository owner:

Name: [Mrinal Arandhara]
Email: [arandharamrinal@gmail.com]

License
This software is distributed under the [INSERT LICENSE HERE]. See the LICENSE file for details.
