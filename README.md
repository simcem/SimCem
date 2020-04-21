# SimCem thermodynamic and process simulation package

This package uses the [stator Comptational Algebra System
(CAS)](https://github.com/toastedcrumpets/stator) to treat
thermodynamics in an equation-oriented way. This is important in
**applied** thermodynamics as often additional constraints such as
heat transfer, kinetics, or multi-stage processes need to be
included.

You may be familiar with Gibb's free energy minimisation, but this on
its own can only consider a single temperature and pressure. Other
constraints need the minimisation of other free energies, for example
including heat loss/gain needs the maximisation of the system entropy.

Using an equation-oriented approach reduces the required iterations to
solve these problems and can lead to better optimisation as the whole
problem is exposed to the optimiser; however, it makes the setting up
of such problems significantly more complex. The idea of SimCem is to
hide as much of this complexity as possible.

Ultimately the thermodynamic problems are reduced to a constrained
minimisation problem which is solved with the
[IPOPT](https://coin-or.github.io/Ipopt) library. Using a CAS allows
these problems to be easily exported into other systems, such as
Matlab/Mathematica, for verification.

# Setup

This library is developed on Linux, but should compile on Windows
(which is still being developed). You will need to install git, Python
3, IPOPT, CMake, and sundials yourself. You can then clone the
repository and its other dependencies (Eigen, sundials, pybind11) like
so:

```
git clone --recurse-submodules http://github.com/toastedcrumpets/SimCem
```

## Current windows issues

The real challenge on windows is getting all the
dependencies. Sundials and IPOPT are the two blockers preventing this
at this time. Once the code is released, CI services such as
Travis/Appveyor will be used to build packages for windows to save
everybody this effort.

# Building

To build the package, run

```
./setup.py build
```

# Installation

To install the package, run

```
./setup.py install
```

# Usage/Database

You will need a thermodynamic database to use SimCem. A free one
containing data from NASA is available on request from
m.campbellbannerman@abdn.ac.uk. This will be released for free in a
separate repository shortly.