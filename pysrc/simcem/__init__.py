from simcem.core import *

# Write some python repr and str implementations
Isotope.__repr__ = lambda s: 'Isotope(symbol='+repr(s.symbol)+', name='+repr(s.name)+', Z='+repr(s.Z)+', N='+repr(s.N)+', mass='+repr(s.mass)+', mass_uncertainty='+repr(s.mass_uncertainty)+', abundance='+repr(s.abundance)+', category='+repr(s.category)+')'
Isotope.__str__ = Isotope.__repr__

Element.__repr__ = lambda s: 'Element(symbol='+repr(s.symbol)+', name='+repr(s.name)+', Z='+repr(s.Z)+', N='+repr(s.N)+', mass='+repr(s.mass)+', mass_uncertainty='+repr(s.mass_uncertainty)+', abundance='+repr(s.abundance)+', category='+repr(s.category)+", group="+repr(s.group)+', period='+repr(s.period)+', block='+repr(s.block)+', referenceComponentID='+repr(s.referenceComponentID)+', isotopes=IsotopeList('+repr([iso for iso in s])+'))'
Element.__str__ = Element.__repr__


#Database.__str__ = lambda s: "Database()"
