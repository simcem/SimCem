# SimCem thermodynamic and process simulation package

[!TIP]
A lot of cement tools based on this library are available online for free here [simcem.streamlit.app](https://simcem.streamlit.app).

This SimCem repository is a thermodynamics package developed by [Marcus bannerman](https://marcusbannerman.co.uk)
and his collaborators in his research on cement production. 

It is a thermodynamic simulation package with a data set suited for calcium sulfoaluminate cement production. Early versions of this software was used in the following papers:
* [I. Garcia, A. Elhoweris, T. Hanein, M. N. Bannerman, and F. P. Glasser, “Advances in clinkering technology of calcium sulfoaluminate cement,” Adv. Cem. Res., 29, 405–417 (2017)](http://dx.doi.org/10.1680/jadcr.17.00028)
* [T. Hanein, I. Galan Garcia, F. P. Glasser, S. Skalamprinos, A. Elhoweris, M.-E. Imbabi, and Bannerman, “Stability of ternesite and the production at scale of ternesite-based clinkers,” Cement Concrete Res., 98, 91-100 (2017) ](http://dx.doi.org/10.1016/j.cemconres.2017.04.010)
* [I. Garcia, T. Hanein, A. Elhoweris, M. N. Bannerman, and F. P. Glasser, “Phase compatibility in the system CaO-SiO2-Al2O3-SO3-Fe2O3 and the effect of partial pressure on phase stability,” Ind. Eng. Chem. Res., 56, (9), 2341-2349 (2017)](http://dx.doi.org/10.1021/acs.iecr.6b03470)
* [T. Hanein, “Development of a novel calcium sulfoaluminate cement production process,” School of Engineering, University of Aberdeen (2016)](https://ethos.bl.uk/OrderDetails.do?uin=uk.bl.ethos.698869)

This library also has a greatly improved implementation of the cement kiln heat transfer simulator listed here:
* [T. Hanein, F. P. Glasser, and M. N. Bannerman, “One-dimensional steady-state thermal model for rotary kilns used in the manufacture of cement,” Adv. Appl. Ceram., 116, (4), 207-215 (2017)](http://dx.doi.org/10.1080/17436753.2017.1303261)

If you're looking for OPC cement production with melt phases and solid solutions, you should use ESPEI and Pycalphad to simulate solid solutions. A TDB database for this is available here:
* [W. Abdul, “Thermodynamic modelling of portland cement clinker production,” School of Engineering, University of Aberdeen (2023)](http://marcusbannerman.co.uk/wahab_abdul_thesis.pdf)
* [W. Abdul, C. Mawalala, A. Pisch, and M. N. Bannerman, “CaO-SiO2 assessment using 3rd generation CALPHAD models,” Cement and Concrete Research, 173, 107309 (2023)](http://dx.doi.org/10.1016/j.cemconres.2023.107309)


## The technical approach 
This package uses the [stator Comptational Algebra System
(CAS)](https://github.com/toastedcrumpets/stator) to treat
thermodynamics in an equation-oriented way. This is important in
applied thermodynamics, as often additional constraints such as heat
transfer, kinetics, or multi-stage processes need to be
included. These are crucial inclusions in process simulation and using
an equation-oriented approach can reduce the required iterations to
solve these types problems and can lead to better optimisation;
however, the setting up of such problems is significantly more
complex. The ambition of SimCem is to hide as much of this complexity
as possible and make process engineering and thermodynamics on real
systems more fun and accessible.

SimCem is still in development so much is aspirational at this
time. The SimCem project also wants to improve thermodynamic fitting
workflows by ensuring data integrity. This means tools are being built
which can take measurements from its sources in original units and
closely linked to the publications/datafiles/graphs to allow careful
and easy review of fitted parameters. This work is mainly the focus of
the website SimCem.com.

## Installation on Linux

Try installing via pip, this should Just Work™ for most Linux/OSX
distributions.

```bash
pip install simcem
``` 

If we've not prebuilt wheels for your system then pip will attempt to
compile simcem but you will need to install cmake, eigen, boost,
sundials, and ipopt for this to succeed. You can get these on Ubuntu
using the following command:

```bash
sudo apt install cmake curl git unzip zip ninja-build g++ gcc make gfortran autoconf libtool python3-dev autoconf-archive pkg-config liblapack-dev libblas-dev
```

## Compiling/Building the package yourself

Make sure to initialise and update submodules 
if you haven't done so with git clone, run this command 

```bash
git clone https://github.com/simcem/simcem
cd simcem
git submodule init
git submodule update
pip install .
```

If you want to run the build yourself, just use CMake directly
```bash
mkdir build
cd build
cmake ..
make
```

## Usage

Please check out the [tests folder](/tests) to see examples of what can be calculated. 
* The file [test_database.py](/tests/test_database.py) shows how to perform combustion calculations and basic phase equilibria.
* [test_kiln_example.py](/tests/test_kiln_example.py) shows how to perform kiln simulations, while [test_kiln.py](/tests/test_kiln.py) iterates over all the kiln data we have from literature and compares our current model's performance.


## Alternatives

SimCem draws inspiration from these projects which are all excellent
solutions and better than SimCem in many ways:

Open/free software
* [Pycalphad](https://pycalphad.org/docs/latest/) is an python Gibbs free energy minimizer which can solve for phase diagrams (among other things). The main advantage over this version of simcem is that it can simulate solid solutions as well as melt phases, so its needed for OPC systems
* [ESPEI](https://espei.org/en/latest/) is built on top of Pycalphad and tries to solve the problems of thermodynamic fitting for metals/melts/solid-solutions systems.
* [OpenCalphad](http://www.opencalphad.com/) a fortran CALPHAD approach
* [BurnMan](https://geodynamics.github.io/burnman/): A nice thermodynamic modelling package for geological processes (so it has pressure-dependence in solids).
* [NASA CEA](https://www1.grc.nasa.gov/research-and-engineering/ceaweb/topicshome/) is a rocketry simulation program and the source of much of the thermodynamic data used in the simcem database. 
