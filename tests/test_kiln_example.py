from simcem import *
from simcem.kiln import *
from scipy.optimize import minimize, leastsq, brenth

# First load the database of properties
db = defaultDatabase()

# Now create a kiln model, this is the KDO at Ibutech, setup for trial point 21
kilnLength = 7.4
kiln = Kiln(RPM=1.8, #Kiln RPM
            innerRadius=0.15, #Inner radius in meters
            length=kilnLength, #Kiln length in meters
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

# Now define the gases
primaryAir = 263 * 1000.0 / 3600.0 #L/s
secondaryAir = 40 * 1000 / 3600 #L/s
methane = 19 * 1000.0 / 3600.0 # L/s

# We now calculate the combustion/flame temperature for the last
# "slice" of the kiln model (at the solids outlet, gas inlet).
# We are calculating the flame temperature of the kiln burner.
init_slice = Slice(solidIn,
                   (primaryAir + secondaryAir) / 1000.0, #Air input in m^3/s
                   methane / 1000.0, #Methane input in m^3/s
                   0 / 1000.0, #O2 input in m^3/s
                   0 / 1000.0, #SO2 input in m^3/s
                   273.15 + 20, #Tsolid K
                   0, #Tgas K (not important, will be calculated anyway)
                   0, #Position of the slice in meters (this is the first slice)
                   db)

# We pull out the temperature of the gas (calculated from an adiabatic
# combustion calculation). This will be our target temperature to get
# the kiln to shoot for.
Tgas_out_target = init_slice.gas.T()


## EXTRA, this is how you'd do everything yourself.
#  #This is an adiabatic flame calculation with a solid fuel
#  inlet_gas = ModelIdealGasTp(db, Components({'CH4':10, 'O2':20, 'N2':80, 'CO':0, 'CO2':0, 'NO':0, 'NO2':0, 'H2O':0}), T=298.15, p=1e5)
#  inlet_solid_fuel = ModelIncompressible(db, Components({"C:graphite":1.0}), 2600, 1.01325e5, "solid", False)
#  sys = System(Objective_t.p, Objective_t.H, True)
#  sys.append(inlet_gas)
#  sys.append(inlet_solid_fuel)
#  sys.equilibrate()
#
## Now we can plug the resulting combusted gas, and the inlet kiln solid into the kiln model.
#  inlet_solid = ModelIncompressible(db, Components({"CaO:Crystal":10}), 298.15, 1.01325e5, "solid", True)
#  init_slice = Slice(inlet_gas, inlet_solid, 0)




# Make a list of points to calculate the kiln conditions at
stop_points = DoubleList()
import numpy as np
#Starting at 0.1, take steps of 0.1 until we reach the kiln length
#Can't start at zero, this would cause an error with the solver.
for v in np.arange(0.1, kilnLength, 0.1): 
    stop_points.append(v)

# Lets make a function that, given a solids inlet temperature of T, returns the difference between our target gas temperature
def func(T):
    #This is a whole solving of a kiln.
    init_slice.gas.set(Objective_t.T, float(T), Objective_t.p)
    kiln.getSlices().clear()
    kiln.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)
    #print("Current exit Tgas=",kiln.getSlices()[len(kiln.getSlices())-1].gas.T())
    #print("Difference is ",kiln.getSlices()[len(kiln.getSlices())-1].gas.T() - T)
    return kiln.getSlices()[-1].gas.T() - Tgas_out_target
#Then solve    
Tgas_in = brenth(func, 298.15, 1200.0, disp=True)

#Here we set the temperature of the initial slice to the solution, and
#run one more calculation to be sure (as the solving above might not
#finish with the "best" solution).
init_slice.gas.set(Objective_t.T, Tgas_in, Objective_t.p)
kiln.getSlices().clear()
kiln.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)


#Grab some data out from the kiln
Z = [slice.Z for slice in kiln.getSlices()] #position
Tg = [slice.gas.T()-273.15 for slice in kiln.getSlices()] #gas temperature
Ts = [slice.solid.T()-273.15 for slice in kiln.getSlices()] #solid temperature
Twall = [slice.T_wall-273.15 for slice in kiln.getSlices()] #wall temperature
Tsh = [slice.T_ext_shell-273.15 for slice in kiln.getSlices()] #Shell temperature

#Lets get the heat fluxes, use a trapezodal rule for integration
lastZ = 0
lastQ = kiln.Q_w_ext(init_slice)
Qtotal = 0
for kslice in kiln.getSlices(): #the slices
    dZ = kslice.Z - lastZ
    newQ = kiln.Q_w_ext(kslice)
    Qtotal += (newQ + lastQ) * dZ / 2
    lastZ = kslice.Z
    lastQ = newQ
print("Total heat loss", Qtotal, 'W')

#We can also check the enthalpy change in the system
Hin = kiln.getSlices()[-1].gas.H() + init_slice.solid.H()
Hout = init_slice.gas.H() + kiln.getSlices()[-1].solid.H()
print("dH =", Hin-Hout, 'W')


#Let's plot the result!
import matplotlib.pyplot as plt
ax = plt.subplot(111)
ax.grid(True)
ax.set_xlim([0.0, kilnLength])
ax.plot(Z, Ts, '-', linewidth=3)
ax.plot(Z, Tg, '-', linewidth=3)
ax.plot(Z, Twall, '--', linewidth=3)
ax.plot(Z, Tsh, '-', linewidth=3)

# Here, we're pulling some experimental data out of the library for
# comparison against
import simcem.KilnData
kiln_data = simcem.KilnData.KDOKilnTP21r()
ax.plot(kiln_data.x_gas_off_wall, [x-273.15 for x in kiln_data.y_gas_off_wall], 'o')
ax.plot(kiln_data.x_bed, [x-273.15 for x in kiln_data.y_bed], 'Dk')
ax.plot(kiln_data.x_gas_off_bed, [x-273.15 for x in kiln_data.y_gas_off_bed], 'pg')
ax.plot(kiln_data.x_wall, [x-273.15 for x in kiln_data.y_wall], 'Hy')

plt.gcf().set_size_inches(5.0, 3.0)
plt.title("Kiln example calculation")
plt.xlabel("Axial Length (m)")
plt.ylabel("Temperature (${}^\\circ$C)")
#plt.show()
