#!/usr/bin/env python3

from collections import defaultdict
import simcem
from simcem import *
import numpy as np

db = simcem.defaultDatabase()

def solve(mode='Taylor'):
    #mass input
    in_mass_comp = simcem.Components({'SiO2':21.5, 'Al2O3':5.2, 'Fe2O3':2.8, 'CaO':66.6, 'MgO':1.0, 'K2O':0.6, 'Na2O':0.2, 'SO3':1.0})
    in_moles_comp = MassToMoles(db, in_mass_comp)
    in_moles_dict = defaultdict(float, in_moles_comp.as_dict())
    
    #Bogue
    if mode == 'Bogue':
        molar_basis = True
        outs = ['Ca3SiO5', 'Ca2SiO4', 'Ca3Al2O6', 'Ca4Fe2Al2O10']
        ins = ['CaO', 'SiO2', 'Al2O3', 'Fe2O3']
        compositions = np.array([
            [3, 1, 0, 0],
            [2, 1, 0, 0],
            [3, 0, 1, 0],
            [4, 0, 1, 1],        
        ])
    elif mode == 'Taylor':
        molar_basis = False
        ins = ['CaO', 'SiO2', 'Al2O3', 'Fe2O3', 'MgO', 'SO3', 'Na2O', 'P2O5', 'K2O', 'TiO2', 'Mn2O3']
        outs = ['Ca3SiO5', 'Ca2SiO4', 'Ca3Al2O6', 'Ca4Fe2Al2O10']
        compositions = np.array([
        [3.000, 0.985, 0.023, 0.01,  0.064, 0,     0.004, 0.003, 0.002, 0,     0,],
        [2.000, 0.926, 0.036, 0.01,  0.022, 0.002, 0.003, 0.002, 0.017, 0.004, 0,],
        [3.000, 0.183, 0.912, 0.095, 0.103, 0,     0.048, 0,     0.022, 0.007, 0,],
        [4.000, 0.283, 1.014, 0.633, 0.351, 0,     0.008, 0,     0.01,  0.095, 0.021,],
        ])
        

    in_moles = np.array([in_moles_dict[comp] for comp in ins])

    #Turn outs (C3S) into oxides (3C, 1S, 0 0)
    #out =  np.dot(np.array([1,0,0,0]), compositions)
    
    out_moles = np.dot(in_moles, np.linalg.pinv(compositions))
    out_mass = [out_moles[idx] * db.getComponent(species).mass() for idx, species in enumerate(outs)]
    print(out_mass)
    print([64.7, 12.9, 9.0, 8.5])
solve()
    
