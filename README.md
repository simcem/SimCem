# SimCem thermodynamic and process simulation package

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

# Installation on Linux/Mac-OSX

Try installing via pip, this should Just Work™ for most Linux/OSX
distributions.

```
pip install simcem
``` 

If we've not prebuilt wheels for your system then pip will attempt to
compile simcem but you will need to install cmake, eigen, boost,
sundials, and ipopt for this to succeed. You can get these on Ubuntu
using the following command:

```
sudo apt-get install libeigen3-dev libsundials-dev coinor-libipopt-dev
```

# Installation on Windows

The only supported method to obtain SimCem is via conda. We package
SimCem occasionally on Windows systems
(here)[https://anaconda.org/toastedcrumpets/simcem]. To install it just run the following:


```
conda install -c toastedcrumpets simcem
```

# Compiling/Building the package yourself

While CMake is used to perform the actual compilation, it is hooked
into python setuptools, so to build the package just run

```
python ./setup.py build
```

# Installation

To install the package, run

```
python ./setup.py install
```

# Usage/Database

You will need a thermodynamic database to use SimCem. A free one
containing data from NASA is available in the [GitHub repository here](https://github.com/toastedcrumpets/SimCem/raw/master/free_database.xml).
Eventually, you will be able to download databases from SimCem.com.

# Alternatives

SimCem draws inspiration from these projects which are all excellent
solutions and better than SimCem in many ways:

Open/free software
* [Pycalphad](https://pycalphad.org/docs/latest/) is an python Gibbs free energy minimizer which can solve for phase diagrams (among other things).
* [ESPEI](https://espei.org/en/latest/) is built on top of Pycalphad and tries to solve the problems of thermodynamic fitting for metals/oslid systems.
* [OpenCalphad](http://www.opencalphad.com/) a fortran CALPHAD approach 
* 
