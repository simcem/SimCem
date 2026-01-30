from simcem.core import *
from simcem.kiln import *
import os

System.setlibHSLpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libhsl.so'))

# Write some python repr and str implementations
Isotope.__repr__ = lambda s: 'Isotope(symbol='+repr(s.symbol)+', name='+repr(s.name)+', Z='+repr(s.Z)+', N='+repr(s.N)+', mass='+repr(s.mass)+', mass_uncertainty='+repr(s.mass_uncertainty)+', abundance='+repr(s.abundance)+', category='+repr(s.category)+')'
Isotope.__str__ = lambda s: 'Isotope(symbol='+repr(s.symbol)+', N='+repr(s.N)+', mass='+repr(s.mass)+')'

Element.__repr__ = lambda s: 'Element(symbol='+repr(s.symbol)+', name='+repr(s.name)+', Z='+repr(s.Z)+', N='+repr(s.N)+', mass='+repr(s.mass)+', mass_uncertainty='+repr(s.mass_uncertainty)+', abundance='+repr(s.abundance)+', category='+repr(s.category)+", group="+repr(s.group)+', period='+repr(s.period)+', block='+repr(s.block)+', referenceComponentID='+repr(s.referenceComponentID)+', isotopes=IsotopeList('+repr([iso for iso in s])+'))'
Element.__str__ = lambda s: 'Element(symbol='+repr(s.symbol)+', mass='+repr(s.mass)+', '+str(len(s))+' isotopes)'

Components.__repr__ = lambda s: 'Components('+repr({k:v for k,v in s.items()})+')'
Components.__str__ = Components.__repr__

####
Database.__str__ = lambda s: 'Database('+str(len(s.getElements()))+' elements, '+str(len(s.getComponents()))+' components)'

Components.as_dict = lambda s: dict(s.items())

def defaultDatabase():
    import os
    return Database(os.path.join(os.path.dirname(__file__), 'free_database.xml'))


def _kilnSolve(self, db, fuelgas, rawmix, fuelSolid=None):
    Tamb=298.15
    sys = System(Objective_t.p, Objective_t.H, True)
    gas = ModelIdealGasTp(db, fuelgas, Tamb, 1.0132e5)
    sys.append(gas)
    if fuelSolid != None:
        sys.append(ModelIncompressible(db, inlet_solid_fuel, Tamb, 1.01325e5, "solid", True))
    sys.equilibrate()

    gasTarget = gas.T()
    
    # Now we can plug the resulting combusted gas, and the inlet kiln solid into the kiln model.
    inlet_solid = ModelIncompressible(db, rawmix, Tamb, 1.01325e5, "solid", True)
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
        init_slice.gas.set(Objective_t.T, float(T), Objective_t.p)
        self.getSlices().clear()
        self.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)
        return self.getSlices()[-1].gas.T() - gasTarget
    #Then solve    
    from scipy.optimize import brenth
    Tgas_in = brenth(func, 298.15, 1200.0, disp=True)
    
    #Here we set the temperature of the initial slice to the solution, and
    #run one more calculation to be sure (as the solving above might not
    #finish with the "best" solution).
    init_slice.gas.set(Objective_t.T, Tgas_in, Objective_t.p)
    self.getSlices().clear()
    self.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)

    return 
    
Kiln.solve = _kilnSolve
