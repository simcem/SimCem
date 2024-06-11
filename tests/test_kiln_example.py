#!/usr/bin/env python3

from simcem import *
from scipy.optimize import minimize, leastsq, brenth
import unittest

# First load the database of properties
def test_kiln_example():
    db = defaultDatabase()

    # Now create a kiln model, this is the KDO at Ibutech, setup for trial point 21
    kiln = Kiln(RPM=1.8, #Kiln RPM
                innerRadius=0.15, #Inner radius in meters
                length=7.4, #Kiln length in meters
                particleDiam=0.0025, #Clinker particle diameter in meters
                shell_emissivity=0.80, #Emissivity of the outer kiln surface
                bed_emissivity=0.9, #Emissivity of the clinker surface
                wall_emissivity=0.85, #Emissivity of the inner kiln wall
                solid_density=2900, #Density of the clinker material (void free) in kg/m^3, this is for CaCO3
                bed_void_frac = 2000 / 2900, #Void fraction of the clinker (bulk density over solid density)
                solid_k=Expr("2.5"), #Thermal conductivity of the clinker, should be a mathematical expression as a function of T
                db=db)

    # We use the simplest bed model, fixed height.
    kiln.setFixedHeightBedModel(0.08) #The loading of the kiln (area/volume fraction filled by clinker)

    # Setup the layers of the kiln wall, giving material name, thickness in meters, and an expression for the thermal conductivity
    kiln.add_layer(material="RefractoryBrick", thickness=0.05, k=Expr("1.39666667 + 2.5e-04 * T"))
    kiln.add_layer(material="SteelShell", thickness=0.005, k=Expr("57"))

    # Lets get started solving the kiln conditions. First, lets define the material inputs

    # Say what solid is coming in, in kg / s (27 kg/hr is given)
    solidIn = Components({'CaCO3':27.0 / 3600.0}) 
    # Then convert it to moles
    solidIn = MassToMoles(db, solidIn)

    # Now define the gases in volumetric flow
    primaryAir = 263 * 1000.0 / 3600.0 #L/s
    secondaryAir = 40 * 1000 / 3600 #L/s
    methane = 19 * 1000.0 / 3600.0 # L/s

    # Convert them into one flow in moles
    ingas = Components({"H2O":0, "CO2":0})
    for flow, gas in [(primaryAir + secondaryAir, Database.getDryAir()), (methane, Components({'CH4':1}))]:
        g = ModelIdealGasTp(db, gas, 298.15, 1.01325e5)
        scale = (flow / g.V()) / 1000
        newgas = ingas + gas * scale
        ingas = newgas

    #Solve the kiln given the fuel and solid inlets!
    kiln.solve(db, ingas, solidIn)

    #Grab some data out from the kiln
    Z = [slice.Z for slice in kiln.getSlices()] #position
    Tg = [slice.gas.T()-273.15 for slice in kiln.getSlices()] #gas temperature
    Ts = [slice.solid.T()-273.15 for slice in kiln.getSlices()] #solid temperature
    Twall = [slice.T_wall-273.15 for slice in kiln.getSlices()] #wall temperature
    Tsh = [slice.T_ext_shell-273.15 for slice in kiln.getSlices()] #Shell temperature

    #Lets get the heat fluxes, use a trapezodal rule for integration
    lastZ = 0
    lastQ = kiln.Q_w_ext(kiln.getSlices()[0])
    Qtotal = 0
    for kslice in kiln.getSlices(): #the slices
        dZ = kslice.Z - lastZ
        newQ = kiln.Q_w_ext(kslice)
        Qtotal += (newQ + lastQ) * dZ / 2
        lastZ = kslice.Z
        lastQ = newQ
    print("Total heat loss", Qtotal, 'W')

    #We can also check the enthalpy change in the system
    Hin = kiln.getSlices()[-1].gas.H() + kiln.getSlices()[0].solid.H()
    Hout = kiln.getSlices()[0].gas.H() + kiln.getSlices()[-1].solid.H()
    print("dH =", Hin-Hout, 'W')


    #Let's plot the result!
    import matplotlib.pyplot as plt
    ax = plt.subplot(111)
    ax.grid(True)
    ax.plot(Z, Ts, '-', linewidth=3)
    ax.plot(Z, Tg, '-', linewidth=3)
    ax.plot(Z, Twall, '--', linewidth=3)
    ax.plot(Z, Tsh, '-', linewidth=3)

    # Here, we're pulling some experimental data out of the library for
    # comparison against
    import simcem.KilnData
    kiln_data = simcem.KilnData.KDOKilnTP21()
    ax.plot(kiln_data.x_gas_off_wall, [x-273.15 for x in kiln_data.y_gas_off_wall], 'o')
    ax.plot(kiln_data.x_bed, [x-273.15 for x in kiln_data.y_bed], 'Dk')
    ax.plot(kiln_data.x_gas_off_bed, [x-273.15 for x in kiln_data.y_gas_off_bed], 'pg')
    ax.plot(kiln_data.x_wall, [x-273.15 for x in kiln_data.y_wall], 'Hy')

    plt.gcf().set_size_inches(5.0, 3.0)
    plt.title("Kiln example calculation")
    plt.xlabel("Axial Length (m)")
    plt.ylabel("Temperature (${}^\\circ$C)")
    #plt.show()

if __name__ == '__main__':
    unittest.main()
