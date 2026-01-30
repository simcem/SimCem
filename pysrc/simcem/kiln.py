import simcem.core as core
from simcem.kiln_core import *

def solve(self, db, fuelgas, rawmix, fuelSolid=None):
    Tamb=298.15
    sys = core.System(core.Objective_t.p, core.Objective_t.H, True)
    gas = core.ModelIdealGasTp(db, fuelgas, Tamb, 1.0132e5)
    sys.append(gas)
    if fuelSolid != None:
        sys.append(core.ModelIncompressible(db, fuelSolid, Tamb, 1.01325e5, "solid", True))
    sys.equilibrate()

    gasTarget = gas.T()
    
    # Now we can plug the resulting combusted gas, and the inlet kiln solid into the kiln model.
    inlet_solid = core.ModelIncompressible(db, rawmix, Tamb, 1.01325e5, "solid", True)
    init_slice = Slice(gas, inlet_solid, 0)
    
    # Make a list of points to calculate the kiln conditions at
    stop_points = DoubleList()
    import numpy as np
    #Starting at 0.1, take steps of 0.1 until we reach the kiln length
    #Can't start at zero, this would cause an error with the solver.
    for v in np.arange(0.1, self.length(), 0.1): 
        stop_points.append(v)

    def func(T):
        #This is a whole solving of a kiln.
        init_slice.gas.set(core.Objective_t.T, float(T), core.Objective_t.p)
        self.getSlices().clear()
        self.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)
        return self.getSlices()[-1].gas.T() - gasTarget
    #Then solve    
    from scipy.optimize import brenth
    Tgas_in = brenth(func, 298.15, 1200.0, disp=True)
    
    #Here we set the temperature of the initial slice to the solution, and
    #run one more calculation to be sure (as the solving above might not
    #finish with the "best" solution).
    init_slice.gas.set(core.Objective_t.T, Tgas_in, core.Objective_t.p)
    self.getSlices().clear()
    self.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)
    return 

Kiln.solve = solve