from simcem import *

db = defaultDatabase()

#### Possible Incompressible phases that can be formed WARNING KEEP OR REMOVE LIQUID PHASES
Incompressiblephases = [
    ('H2O', "water"),
    ('H2SO4', "Liquid"),
    ('NaAlSi3O8', "Liquid"),
    ('NaAlSi3O8', "albite"),
    ('NaAlSi3O8', "albite (high)"),
    ('NaFeSi2O6', "acmite"),
    ('NaAlSi2O6', "jadeite"),
    ('CaAl2Si2O8', "Liquid"),
    ('CaMgSi2O6', "Liquid"),
    ('CaMgSi2O6', "diopside"),
    ('CaMg(CO3)2', "dolomite"),
    ('Mg2Si2O6', "Liquid"),
    ('Mg2Si2O6', "enstatite"),
    ('Fe2SiO4', "Liquid"),
    ('Fe2SiO4', "fayalite"),
    ('Fe2Si2O6', "ferrosilite"),
    #('K3Fe0.5Al4Si19.5O47', "Liquid"),
    #('K3Mg0.5Al4Si19.5O47', "Liquid"),
    ('KAlSi3O8', "Liquid"),
    ('KAlSi3O8', "sanidine"),
    ('Al2SiO5', "Liquid"),
    ('Fe2O3', "hematite"),
    ('SiO2', "alpha-quartz"), 
    ('SiO2', "beta-quartz"),
    ('MgCO3', "magnesite"),
    ('KAlO2', "II"),
    ('KAlO2', "I"),
    ('Na2CO3', "a"),
    ('Na2CO3', "b"),
    ('Na2CO3', "c"),
    ('Na2CO3', "Liquid"),
    ('NaAlO2', "a"),
    ('NaAlO2', "b"),
    ('K2O', "c"),
    ('K2O', "b"),
    ('K2O', "a"),
    ('K2O', "Liquid"),
    ('Na2O', "c"),
    ('Na2O', "b"),
    ('Na2O', "a"),
    ('Na2O', "Liquid"),
    ('Na2SO4', "V"),
    ('Na2SO4', "IV"),
    ('Na2SO4', "I"),
    ('Na2SO4', "Liquid"),
    ('K2CO3', "a"),
    ('K2CO3', "b"),
    ('K2CO3', "Liquid"),
    ('K2SO4', "II"),
    ('K2SO4', "I"),
    ('K2SO4', "Liquid"),
    ('K2Si2O5', "a"),
    ('K2Si2O5', "b"),
    ('K2Si2O5', "c"),
    ('K2Si2O5', "Liquid"),
    ('K2SiO3', "Crystal"),
    ('K2SiO3', "Liquid"),  
    ('MgCO3', "Liquid"),
    ('MgSiO3', "I"),
    ('MgSiO3', "II"),
    ('MgSiO3', "III"),
    ('MgSiO3', "Liquid"),
    ('MgFe2O4', "alpha"),
    ('MgFe2O4', "beta"),
    ('MgFe2O4', "gamma"),
    ('MgAl2O4', "Crystal"),
    ('Mg2SiO4', "Crystal"),
    ('MgO', "Crystal"),
    ('CaAl2SiO6', "Ca-Al clinopyroxene"),
    ('CaFe2O4', "Liquid"),
    ('CaFe2O4', "Crystal"),
    ('Ca4Fe2Al2O10', "Crystal"),
    ('Al2O3', "alpha"),
    ('Ca3Al2Si3O12', "grossular"),
    ('CaSiO3', "cyclowollastonite"),
    ('CaSiO3', "wollastonite"), 
    ('Al2SiO5', "sillimanite"),
    ('Al2SiO5', "andalusite"),
    ('Al2SiO5', "kyanite"),
    ('CaAl2Si2O8', "anorthite"),
    ('CaAl4Si2O10(OH)2', "margarite"),
    ('Ca2Al3Si3O12(OH)', "zoisite"),
    ('Ca2Al2Si3O10(OH)2', "prehnite"),
    ('Ca2SiO4', "beta"),
    ('Al2Si4O10(OH)2', "pyrophyllite"),
    ('HAlO2', "boehmite"),
    ('HAlO2', "diaspore"),
    ('Ca2Al2SiO7', "gehlenite"),
    ('CaAl12O19', "Crystal"),
    #('Al2TiO5', "Crystal"), From Barin Before
    #('Fe2TiO4', "Crystal"),From Barin Before
    #('Fe2TiO5', "Crystal"),From Barin Before
    ('MgTiO3', "Crystal"),
    ('Mg2TiO4', "Crystal"),
    ('MgTi2O5', "Crystal"),
    ('Ti3O5', "a"),
    ('Ti3O5', "b"),
    ('TiO2', "Crystal"),
    ('TiO', "a"),
    ('TiO', "b"),
    ('TiO', "c"),
    ('Ti2O3', "I"),
    ('Ti2O3', "I'"),
    ('Ti4O7', "Crystal"),
    #('Zn2TiO4', "Crystal"),
    ('CaTiO3', "alpha"),
    ('CaTiO3', "beta"),
    ('FeTiO3', "Crystal"),
    ('FeTiO3', "Liquid"),
    ('Ca2Fe2O5', "Liquid"), 
    ('Ca2Fe2O5', "Crystal"),
    ('Al6Si2O13', "mullite"),
    ('Ca5Si2C2O13', "tilleyite"),   
    ('Ca12Al14O33', "alpha"),
    ('Ca12Al14O33', "beta"),
    ('CaAl2O4', "krotite"),
    ('CaAl4O7', "Crystal"),  
    ('Ca3Al2O6', "Crystal"),
    ('Ca5Si2CO11', "spurrite"),
    ('Ca2SiO4', "alphaprime"),
    ('Ca2SiO4', "gamma"),
    ('Ca2SiO4', "beta"),
    ('Ca2SiO4', "alpha"),
    ('Ca3SiO5', "Crystal"),
    ('Ca5Si2SO12', "Crystal"),
    ('Ca4Al6SO16', "pseudocubic"),
    ('Ca4Al6SO16', "orthorhombic"),
    #('Ca4Al5.5Fe0.5SO16', "Yeelimite_F"),
    ('Ca3Si2O7', "rankinite"),                ########## CHECK THIS !!!!!!! WHICH RANKINITE IS IT?
    ('CaO', "Crystal"),
    ('CaCO3', "Crystal"),
    ('CaCO3', "Liquid"), 
    ('CaSO4', "II"),
    ('CaSO4', "I"),
    ('CaSO4', "Liquid"),
    ]

solids=[]
liquids=[]

for species,phase in Incompressiblephases:
    ID = species+":"+phase
    phaseref = db.getComponent(species).getPhase(phase, '')
    if not len(phaseref):
        print("No data available for '"+species+"' '"+phase)
        
    if phaseref.type == 'solid':
        solids.append((species, phase))
    elif phaseref.type == 'liquid':
        liquids.append((species, phase))
    else:
        raise Exception("Failed to find a matching solid or liquid phase for '"+species+"' '"+phase+"'\nFound phases"+str(db.getComponent(species)))

def setup_phases(gas_components, in_components, T):
    #We try inserting components as solids and liquids
    to_insert = dict(in_components.items()) # Force a copy to prevent changing the original in_components

    solid_comps = {}
    for species, phase in solids:
        ID=species+':'+phase
        testphase = ModelIncompressible(db, Components({ID:1}), T, 1.01325e5, "solid", False)
        if species in to_insert:
            if testphase.hasData():
                solid_comps[ID] = to_insert[species]
                del to_insert[species]
        else:
            solid_comps[ID] = 0
            
    liquid_comps = {}
    for species, phase in liquids:
        ID=species+':'+phase
        testphase = ModelIncompressible(db, Components({ID:1}), T, 1.01325e5, "liquid", False)
        if species in to_insert:
            if testphase.hasData():
                liquid_comps[ID] = to_insert[species]
                del to_insert[species]
        else:
            liquid_comps[ID] = 0

    if to_insert:
        raise Exception("Failed to insert solid species as solid or liquid,"+str(to_insert))

    #Create the system
    sys = System(Objective_t.p, Objective_t.T, True, max_eval=1000)

    #Now create the phases
    solid  = ModelIncompressible(db, Components(solid_comps), T, 1.01325e5, "solid", mixing=False)
    liquid = ModelIncompressible(db, Components(liquid_comps), T, 1.01325e5, "liquid", mixing=False)
    gas    = ModelIdealGasTp(db, gas_components, T, p=1.01325e5)
    
    if len(gas_components) > 0:
        sys.append(gas)

    if len(liquid_comps) > 0:
        sys.append(liquid)
    
    sys.append(solid)
    return gas, solid, liquid, sys

def clinkerize(gas_components, in_components, T):
    gas, solid, liquid, sys = setup_phases(gas_components, in_components, T)

    try:
        sys.equilibrate()
    except Exception as e: 
        print("Did not converge for", gas_components, in_components, T)
    return gas, solid, liquid, sys


def boguey_calcs(input, ins, outs, compositions_inv):
    #Convert mass into moles
    input_components_mass = Components(input) # 'P2O5':P2O5, 'Mn2O3':Mn2O3, 'SO3':SO3, , 'Na2O':Na2O, 'TiO2':TiO2 
    input_components_moles = MassToMoles(db, input_components_mass)
    #Make a lookup that defaults to zero for missing components
    input_moles_dict = defaultdict(float, input_components_moles.as_dict())
    #Now create the vector of inputs
    in_moles = np.array([input_moles_dict[comp] for comp in ins])
    #Create the vector of outputs
    out_moles = np.dot(in_moles, compositions_inv)
    #Make the output dictionary
    output_components_mass = defaultdict(float, {key:val * db.getComponent(key).mass() for key,val in zip(outs, out_moles)})
    return output_components_mass


import numpy as np
from collections import defaultdict
import math

def bogue_calculation(input):
    bogue_outs = ['Ca3SiO5', 'Ca2SiO4', 'Ca3Al2O6', 'Ca4Fe2Al2O10']
    bogue_ins = ['CaO', 'SiO2', 'Al2O3', 'Fe2O3']
    bogue_compositions = np.array([
                [3, 1, 0, 0],
                [2, 1, 0, 0],
                [3, 0, 1, 0],
                [4, 0, 1, 1],        
    ])
    bogue_compositions_inv = np.linalg.pinv(bogue_compositions)
    return boguey_calcs(input, bogue_ins, bogue_outs, bogue_compositions_inv)

def taylor_calculation(input):
    taylor_ins = ['CaO', 'SiO2', 'Al2O3', 'Fe2O3', 'MgO', 'SO3', 'Na2O', 'P2O5', 'K2O', 'TiO2', 'Mn2O3']
    taylor_outs = ['Ca3SiO5', 'Ca2SiO4', 'Ca3Al2O6', 'Ca4Fe2Al2O10']
    taylor_compositions = np.array([
        [3.000, 0.985, 0.023, 0.01,  0.064, 0,     0.004, 0.003, 0.002, 0,     0,],
        [2.000, 0.926, 0.036, 0.01,  0.022, 0.002, 0.003, 0.002, 0.017, 0.004, 0,],
        [3.000, 0.183, 0.912, 0.095, 0.103, 0,     0.048, 0,     0.022, 0.007, 0,],
        [4.000, 0.283, 1.014, 0.633, 0.351, 0,     0.008, 0,     0.01,  0.095, 0.021,],
    ])
    taylor_compositions_inv = np.linalg.pinv(taylor_compositions)
    return boguey_calcs(input, taylor_ins, taylor_outs, taylor_compositions_inv)
    
if True:
    # Unit test of the bogue equations. From https://www.understanding-cement.com/bogue.html
    # Note the subtraction of 1% free lime
    input = {'SiO2':21.5, 'Al2O3':5.2, 'Fe2O3': 2.8, 'CaO':66.6-1.0, 'MgO':1.0, 'K2O':0.6, 'Na2O':0.2, 'SO3':1.0}
    expected_output = {'Ca3SiO5':64.7, 'Ca2SiO4':12.9, 'Ca3Al2O6':9.0, 'Ca4Fe2Al2O10':8.5}
    result = bogue_calculation(input)
    for key, value in expected_output.items():
        assert(math.isclose(result[key], expected_output[key], abs_tol=0.2))