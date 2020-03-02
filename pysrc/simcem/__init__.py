from simcem.core import *

# Write some python representations
Isotope.__repr__ = lambda s: 'Isotope(Z='+repr(s.Z)+', N='+repr(s.N)+', mass='+repr(s.mass)+', mass_uncertainty='+repr(s.mass_uncertainty)+', abundance='+repr(s.abundance)+', symbol='+repr(s.symbol)+', name='+repr(s.name)+', category='+repr(s.category)+')'
Isotope.__str__ = Isotope.__repr__

Database.__str__ = lambda s: "Database()"
