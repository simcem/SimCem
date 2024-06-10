#!/usr/bin/env python3
import os, sys, numpy
import matplotlib.pyplot as plt
import simcem 

db = simcem.defaultDatabase()

F=96485.3329

def addPhases(sys, comps, T, p, singlephase, phaset):
    '''A convenience function to add components to a system, either as a set of
    individual phases (mix=False) or as a single phase
    (mix=True). phaset is a functor that makes the kind of phase you
    want.
    '''
    if len(comps) != 0:
        if singlephase:
            sys.append(phaset(db, simcem.Components(comps), T, p))
        else:
            for c,v in comps.items():
                sys.append(phaset(db, simcem.Components({c:v}), T, p))

def getSystem(solids = {}, gases = {}, p = 1.013e5, T=298.15, mixgases=True):
    '''Given a dictionary of solids/gases with their amounts, it creates
    a thermodynamic system containing them using the NASA model.'''
    sys = simcem.System(simcem.Objective_t.p, simcem.Objective_t.T, False)
    addPhases(sys, gases,  T, p, singlephase=mixgases, phaset=simcem.ModelIdealGasTp)
    addPhases(sys, solids, T, p, singlephase=True, phaset=lambda a,b,c,d : simcem.ModelIncompressible(a,b,c,d, "solid", False))
    return sys
    
def deltaG(solidIn, gasIn, solidOut, gasOut, T, mixgases=False, p=1.013e5):
    '''Calculate a deltaG for a set of gases and solids.'''
    Reactants = getSystem(solidIn,  gasIn,  p=p, T=T, mixgases=mixgases)
    Products  = getSystem(solidOut, gasOut, p=p, T=T, mixgases=mixgases)
    return Products.G() - Reactants.G()

##  sys = getSystem(solids={"CaCO3:aragonite": 1}, gases={"CO2":1},  p=1.013e5, T=298.15, mixgases=True)
##  print("#### ", sys.A(), sys.H()/sys.N(), sys.Cp()/sys.M())
##  
##  ## Here we're just creating plots of deltaG and dividing it by
##  ## Faradays constant (F), and the number of electrons involved, to
##  ## give a voltage
##  T1 = 900+273.15
##  T2 = 1000+273.15
##  p=1.013e5
##  Ts=numpy.arange(0+273.15, 4000+273.15, 100)
##  reac2 = [deltaG({}, {"CO2":1}, {}, {"CO":1, "O2":0.5}, T)/F/2 for T in Ts]
##  reac3 = [deltaG({}, {"H2O":1}, {}, {"H2":1, "O2":0.5}, T)/F/2 for T in Ts]
##  Ts = [T-273.15 for T in Ts]
##  plt.plot(Ts, reac2, '.-k', label=r"$\mathrm{CO}_2\,\to\,\mathrm{CO}+\frac{1}{2}\mathrm{O}_2$")
##  plt.plot(Ts, reac3, '.-r', label=r"$\mathrm{H}_2\mathrm{O}\,\to\,\mathrm{H}_2+\frac{1}{2}\mathrm{O}_2$")
##  plt.legend()
##  plt.ylabel("Potential [V]")
##  plt.ylim(0, 1.3)
##  plt.xlabel("Temperature [C]")
##  plt.show()
##  
##  Ts=numpy.arange(0+273.15, 1300+273.15, 100)
##  reac1 = [deltaG({}, {"CO2":1}, {}, {"CO":1, "O2":0.5}, T)/F/2 for T in Ts]
##  reac2 = [deltaG({"Li2CO3":1}, {}, {"Li2O":1}, {"CO":1, "O2":0.5}, T)/F/2 for T in Ts]
##  reac3 = [deltaG({"Li2CO3":1}, {}, {"Li2O":1, "C:graphite":1.0}, {"O2":1}, T)/F/4 for T in Ts]
##  #No idea what the electrochemistry is of this, guessing 4?
##  reac4 = [deltaG({"Li2O":1}, {"CO2":1}, {"Li2CO3":1}, {}, T)/F/4 for T in Ts]
##  
##  Ts = [T-273.15 for T in Ts]
##  plt.plot(Ts, reac1, '.-r',  label=r"$\mathrm{CO}_2\,\to\,\mathrm{CO}+\frac{1}{2}\mathrm{O}_2$")
##  plt.plot(Ts, reac2, 's-b', label=r"$\mathrm{Li}_2\mathrm{CO}_3\,\to\,\mathrm{CO}+\mathrm{Li}_2\mathrm{O}+\frac{1}{2}\mathrm{O}_2$")
##  plt.plot(Ts, reac3, '.-b', label=r"$\mathrm{Li}_2\mathrm{CO}_3\,\to\,\mathrm{C}+\mathrm{Li}_2\mathrm{O}+\mathrm{O}_2$")
##  plt.plot(Ts, reac4, '.-g', label=r"$\mathrm{Li}_2\mathrm{O}+\mathrm{CO}_2\,\to\,\mathrm{Li}_2\mathrm{CO}_3$")
##  plt.legend()
##  plt.ylabel("Potential [V]")
##  plt.xlabel("Temperature [C]")
##  #plt.ylim(0.7, 2.3)
##  plt.show()
##  
##  T=25+273.15
##  print("E_CO2 split (with mix) = ", deltaG({}, {"CO2":1}, {}, {"CO":1, "O2":0.5}, T, mixgases=True)/F/2, "V")
##  print("E_CO2 split (w/o mix)  = ", deltaG({}, {"CO2":1}, {}, {"CO":1, "O2":0.5}, T, mixgases=False)/F/2, "V")


#### Assumptions with OxyFuel
# We assume that the melt is in equilibrium with the capture gas, and
# that products are removed efficiently from around the cathode/anode.

Ts=numpy.arange(0+273.15, 1300+273.15, 100)

def equilAtmosCalc(solidsIn, solidsOut, gasOut, T, p, gasInCompMoles):
    #Determine the equilibrium gas in
    sys = simcem.System(simcem.Objective_t.p, simcem.Objective_t.T, True, max_eval=1000)
    gasIn = simcem.ModelIdealGasTp(db, gasInCompMoles, T, p)
    sys.append(gasIn)
    sys.equilibrate()
    
    Reactants = simcem.System(simcem.Objective_t.p, simcem.Objective_t.T, False)
    Reactants.append(simcem.ModelIdealGasTp(db, gasIn.components, T, p))
    Reactants.append(simcem.ModelIncompressible(db, simcem.Components(solidsIn), T, p, "solid", False))

    Products = simcem.System(simcem.Objective_t.p, simcem.Objective_t.T, True, max_eval=1000)
    gasOutModel = simcem.ModelIdealGasTp(db, gasIn.components + simcem.Components(gasOut), T, p)
    Products.append(gasOutModel)
    Products.equilibrate() #Equilibrate the gas phase separately??
    Products.append(simcem.ModelIncompressible(db, simcem.Components(solidsOut), T, p, "solid", False))
    
    print("Exit Gas", gasOutModel)
    return Products.G() - Reactants.G()

def oxyfuelCalc(solidsIn, solidsOut, gasOut,p=1.013e5):
    #Oxyfuel capture stream
    gasInCompMass = simcem.Components({'CO2':16.3599, 'O2':0.8064, 'N2':3.4971, "CO":0})
    gasInCompMoles = 10000*simcem.MassToMoles(db, gasInCompMass)
    data = []
    for T in Ts:
        try:
            data.append(equilAtmosCalc(solidsIn, solidsOut, gasOut, T, p, gasInCompMoles))
        except:
            data.append(numpy.nan)
    return numpy.array(data)

#Methane combustion at 2600 with excess of O2 from our test examples
sys = simcem.System(simcem.Objective_t.p, simcem.Objective_t.T, True, max_eval=1000)
combustionGas = simcem.ModelIdealGasTp(db, simcem.Components({"O2":0.21, "N2":0.79, "CO":0, "CO2":0, "H2O":0, "CH4":0.1}), T=273.15+600, p=1.013e5)
sys.append(combustionGas)
sys.equilibrate()
print("COMBUSTION GAS", combustionGas)
#print("COMBUSTION GAS", combustionGas.components / combustionGas.N())
def postCombustion(solidsIn, solidsOut, gasOut, p=1.013e5):
    gasInCompMoles = 10000 * combustionGas.components
    #simcem.Components({"CH4":8.688777086942181e-12, "CO":0.03453317988062815, "CO2":0.07046682011068307, "N2":0.79, "O2":0.01726658995769161}) #"H2O":0.2099999999826224,
    data = []
    for T in Ts:
        try:
            data.append(equilAtmosCalc(solidsIn, solidsOut, gasOut, T, p, gasInCompMoles))
        except:
            data.append(numpy.nan)
    return numpy.array(data)

def atmospheric(solidsIn, solidsOut, gasOut, p=1.013e5):
    #Methane combustion at 2600 with excess of O2 from our test examples
    gasInCompMoles = 10000 * simcem.Components({"N2":0.79, "CO":0.0, "CO2":400e-6, "O2":0.21-0.000500}) #"H2O":0.2099999999826224,
    print("Atmospheric DATA")
    data = []
    for T in Ts:
        try:
            data.append(equilAtmosCalc(solidsIn, solidsOut, gasOut, T, p, gasInCompMoles))
        except:
            data.append(numpy.nan)
    return numpy.array(data)

gasInCompMass = simcem.Components({'CO2':16.3599, 'O2':0.8064, 'N2':3.4971, "CO":0})
gasInCompMoles = simcem.MassToMoles(db, gasInCompMass)
print(gasInCompMoles * 100 / gasInCompMoles.N())

gasInCompMoles = simcem.Components({"CH4":8.688777086942181e-12, "CO":0.03453317988062815, "CO2":0.07046682011068307, "N2":0.79, "O2":0.01726658995769161}) #"H2O":0.2099999999826224,
print(100 * gasInCompMoles / gasInCompMoles.N())


for title, state in [("Oxyfuel post-combustion atmosphere", oxyfuelCalc), ("Post-combustion atmosphere", postCombustion), ("Experiment", atmospheric)]:
    reac1 = state({}, {}, {"CO2":-1, "O2":0.5, "CO":1})/F/2
    reac2 = state({"Li2CO3":1}, {"Li2O":1}, {"CO":1, "O2":0.5})/F/2
    reac3 = state({"Li2CO3":1}, {"Li2O":1, "C:graphite":1.0}, {})/F/4
    #No idea what the electrochemistry is of this, guessing 4?
    reac4 = state({"Li2O":1}, {"Li2CO3":1}, {"CO2":-1})/F/4
    reac5 = state({"Li2CO3":2}, {"Li2O":1, "Li2C2":1}, {"O2":2}) /F/10

    Tplot = [T-273.15 for T in Ts]
    plt.plot(Tplot, reac1, '.-r',  label=r"$\mathrm{CO}_2\,\to\,\mathrm{CO}+\frac{1}{2}\mathrm{O}_2$")
    plt.plot(Tplot, reac2, 's-b', label=r"$\mathrm{Li}_2\mathrm{CO}_3\,\to\,\mathrm{CO}+\mathrm{Li}_2\mathrm{O}+\frac{1}{2}\mathrm{O}_2$")
    plt.plot(Tplot, reac3, '.-b', label=r"$\mathrm{Li}_2\mathrm{CO}_3\,\to\,\mathrm{C}+\mathrm{Li}_2\mathrm{O}+\mathrm{O}_2$")
    plt.plot(Tplot, reac4, '.-g', label=r"$\mathrm{Li}_2\mathrm{O}+\mathrm{CO}_2\,\to\,\mathrm{Li}_2\mathrm{CO}_3$")
    plt.plot(Tplot, reac5, '.-k', label=r"$2\,\mathrm{Li}_2\mathrm{CO}_3 \to \mathrm{Li}_2\mathrm{C}_2 + \mathrm{Li}_2\mathrm{O}+2\,\mathrm{O}_2$")
    plt.legend()
    plt.ylabel("Potential [V]")
    plt.xlabel("Temperature [°C]")
    plt.title(title)
    plt.ylim(-0.5, 2.5)
    plt.show()
