#!/usr/bin/env python3

import simcem
from simcem.clinker import clinkerize
from scipy.optimize import minimize, leastsq, brenth, LinearConstraint
import numpy as np

db = simcem.defaultDatabase()

T= 1125
C4AF = simcem.ModelIncompressible(db, simcem.Components({"Ca4Fe2Al2O10:Crystal": 1}), T, 1.01325e5, "solid", mixing=False)
C = simcem.ModelIncompressible(db, simcem.Components({"CaO:Crystal": 1}), T, 1.01325e5, "solid", mixing=False)
A = simcem.ModelIncompressible(db, simcem.Components({"Al2O3:alpha": 1}), T, 1.01325e5, "solid", mixing=False)
F = simcem.ModelIncompressible(db, simcem.Components({"Fe2O3:hematite": 1}), T, 1.01325e5, "solid", mixing=False)

C2F = simcem.ModelIncompressible(db, simcem.Components({"Ca2Fe2O5:Crystal": 1}), T, 1.01325e5, "solid", mixing=False)

print(C4AF.H()-4*C.H()-A.H()-F.H())
print(C2F.H() - 2*C.H() - F.H(), C2F.S())
