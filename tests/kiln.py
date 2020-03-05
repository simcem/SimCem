import unittest
import time

from simcem import *
from simcem.kiln import *
from scipy.optimize import minimize, leastsq, brenth
import simcem.KilnData
import pylab as plt

class CoreTest(unittest.TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3fms" % (self.id(), t*1000))
    
    @classmethod
    def setUpClass(cls):
        cls.db = Database("free_database.xml")

    def test_kiln(self):
        kiln_exp_data = simcem.KilnData.getExperimentalKilns()
        for kiln_data in kiln_exp_data:
            kiln_data.finish()
            
            kiln = Kiln(RPM=kiln_data.RPM,
                        innerRadius=kiln_data.innerRadius,
                        length=kiln_data.length,
                        particleDiam=kiln_data.particleDiameter,
                        shell_emissivity=kiln_data.shellEmissivity,
                        bed_emissivity=kiln_data.bedEmissivity,
                        wall_emissivity=kiln_data.wallEmissivity,
                        solid_density=kiln_data.solid_density,
                        bed_void_frac = kiln_data.bulk_density / kiln_data.solid_density,
                        solid_k=kiln_data.k_expr,
                        db=self.db)
            kiln.setFixedHeightBedModel(kiln_data.solidLoading)
            
            for layer in kiln_data.insulationLayers:
                kiln.add_layer(material=layer[0], thickness=layer[1], k=Expr(layer[2]))    
        
            #Gas off bed appears to have large fluctuations, probably due to
            #issues with keeping it clear of the bed.
            #+kiln_data.x_gas_off_bed
            stop_points_py = list(set(kiln_data.x_gas_off_wall+kiln_data.x_bed+kiln_data.x_wall+[kiln_data.length]))
            stop_points_py.sort()
            stop_points = DoubleList()
            for v in stop_points_py:
                stop_points.append(v)
        
            if hasattr(kiln_data, 'SolveMode'):
                init_slice = Slice(MassToMoles(self.db, kiln_data.solidIn),
                                   (kiln_data.primaryAir + kiln_data.secondaryAir) / 1000.0, #m^3/s
                                   kiln_data.inGasFlow / 1000.0, #m^3/s
                                   kiln_data.inO2Flow / 1000.0, #m^3/s
                                   kiln_data.inSO2Flow / 1000.0, #m^3/s
                                   273.15 + 20, #Tsolid K
                                   0, #Tgas K
                                   0, #Z0
                                   self.db)
                Tgas_out_target = init_slice.gas.T()
                #print("Solving for T_gas_out=",init_slice.gas.T())
                def func(T):
                    #print("Trying Tgas_in=", T)
                    init_slice.gas.set(Objective_t.T, float(max(T,298.15)), Objective_t.p)
                    kiln.getSlices().clear()
                    kiln.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)
                    #print("Current exit Tgas=",kiln.getSlices()[len(kiln.getSlices())-1].gas.T())
                    #print("Difference is ",kiln.getSlices()[len(kiln.getSlices())-1].gas.T() - T)
                    return kiln.getSlices()[len(kiln.getSlices())-1].gas.T() - Tgas_out_target
        
                Tgas_in = brenth(func, 298.15, 1200.0, disp=True)
                init_slice.gas.set(Objective_t.T, Tgas_in, Objective_t.p)
                kiln.getSlices().clear()
                kiln.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)
            else:
                try:
                    init_slice = Slice(MassToMoles(self.db, kiln_data.solidIn),
                                       (kiln_data.primaryAir + kiln_data.secondaryAir) / 1000.0, #m^3/s
        		               kiln_data.inGasFlow / 1000.0, #m^3/s
        		               0, #m^3/s
        		               0, #m^3/s
        		               Tsolid = 20+273.15, #Tsolid K
        		               Tgas = 800, #Tgas K
                                       Z0 = 0,
                                       db=self.db)
                except Exception as e:
                    print("Failed to set up initial slice", kiln_data.kilnName+kiln_data.ID)
                    continue
                
                def func(T):
                    init_slice.gas.set(Objective_t.T, float(max(T[0],298.15)), Objective_t.p)
                    init_slice.solid.set(Objective_t.T, float(max(T[1],298.15)), Objective_t.p)
                    #print(init_slice.gas, init_slice.solid)
                    kiln.getSlices().clear()
                    try:
                        kiln.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)
                    except Exception as e:
                        print(kiln_data.kilnName+kiln_data.ID, " FAILED, trying again")
                        kiln.solve_SS_inert(init_slice=init_slice, stop_points=stop_points, store_intermediate=True)
                        return 
                    results = {}
                    for slice in kiln.getSlices():
                        results[slice.Z] = [slice.gas.T(), slice.solid.T(), slice.T_wall]
                
                    error =[]
                    #error = 0
                    #See above for why it has been removed
                    #for Z,T in zip(kiln_data.x_gas_off_bed, kiln_data.y_gas_off_bed):
                    #    error += (results[Z][0] - T)*(results[Z][0] - T)
                    
                    for Z,T in zip(kiln_data.x_gas_off_wall, kiln_data.y_gas_off_wall):
                        error.append(results[Z][0] - T)
                        #error += (results[Z][0] - T)*(results[Z][0] - T)
                    
                    for Z,T in zip(kiln_data.x_bed, kiln_data.y_bed):
                        error.append(results[Z][1] - T)
                        #error += (results[Z][1] - T)*(results[Z][1] - T)
                
                    for Z,T in zip(kiln_data.x_wall, kiln_data.y_wall):
                        error.append(results[Z][2] - T)
                        #error += (results[Z][2] - T)*(results[Z][2] - T)
                        
                    return error
                
                #Sol = minimise(func, [min(kiln_data.y_gas_off_wall),min(kiln_data.y_bed)],method="nelder-mead",options={'ftol': 1,'xtol': 0.2,'disp': True}, bounds=[(298.15, None), (298.15, None)])
                Sol = leastsq(func, [min(kiln_data.y_gas_off_wall),min(kiln_data.y_bed)],ftol=0.001,xtol=0.01)
                func(Sol[0])
            
            #for slice in kiln.getSlices():
            #    kiln.printFluxes(slice)
            
            Z = [slice.Z for slice in kiln.getSlices()]
            Tg = [slice.gas.T()-273.15 for slice in kiln.getSlices()]
            Ts = [slice.solid.T()-273.15 for slice in kiln.getSlices()]
            Twall = [slice.T_wall-273.15 for slice in kiln.getSlices()]
            Tsh = [slice.T_ext_shell-273.15 for slice in kiln.getSlices()]
            
            plt.clf()
            ax = plt.subplot(111)
            ax.grid(True)
            ax.set_xlim([0.0,kiln_data.length])
            ax.plot(Z, Ts, '-', linewidth=3)
            ax.plot(Z, Tg, '-', linewidth=3)
            ax.plot(Z, Twall, '--', linewidth=3)
            ax.plot(Z, Tsh, '-', linewidth=3)
            ax.plot(kiln_data.x_gas_off_wall, [x-273.15 for x in kiln_data.y_gas_off_wall], 'o')
            ax.plot(kiln_data.x_bed, [x-273.15 for x in kiln_data.y_bed], 'Dk')
            ax.plot(kiln_data.x_gas_off_bed, [x-273.15 for x in kiln_data.y_gas_off_bed], 'pg')
            ax.plot(kiln_data.x_wall, [x-273.15 for x in kiln_data.y_wall], 'Hy')
        
            box = ax.get_position()
            #ax.set_position([box.x0, box.y0, box.width * 0.83, box.height])
            factor=1.0
            plt.gcf().set_size_inches(5.0*factor, 3.0*factor)
        
            #ax.legend(["$T_s$", "$T_g$", "$T_w$", "$T_{sh}$","$T^*_{g-off-wall}$","$T^*_{s}$","$T^*_{g-off-bed}$" ,"$T^*_{w}$" ],loc='center left', bbox_to_anchor=(1, 0.5),fancybox=True, shadow=True)
        
            plt.title("Experiment "+kiln_data.kilnName+kiln_data.ID)
            plt.xlabel("Axial Length (m)")
            plt.ylabel("Temperature (${}^\\circ$C)")
            plt.savefig(kiln_data.kilnName+kiln_data.ID+".pdf", bbox_inches="tight", pad_inches=0, transparent=True)
            print(kiln_data.kilnName+kiln_data.ID, " Done")

        
#Fitting of the Tscheng Kiln insulation
#
#
#
#kiln_data = KilnData.TschengKilnA16()
#kiln = Kiln(RPM=kiln_data.RPM,
#            innerRadius=kiln_data.innerRadius,
#            length=kiln_data.length,
#            particleDiam=kiln_data.particleDiameter,
#            shell_emissivity=kiln_data.shellEmissivity,
#            bed_emissivity=kiln_data.bedEmissivity,
#            wall_emissivity=kiln_data.wallEmissivity,
#            solid_density=2627.0,#TScheng value for ottawa sand
#            solid_k=Expr("0.268"),#TScheng value! 
#            db=db)
#kiln.setFixedHeightBedModel(kiln_data.solidLoading)
#for layer in kiln_data.insulationLayers:
#    kiln.add_layer(material=layer[0], thickness=layer[1], k=Expr(layer[2]))    
#slice = Slice(solid=MassToMoles(db, kiln_data.solidIn),
#              volAirFlow = (kiln_data.primaryAir + kiln_data.secondaryAir) / 1000.0, #m^3/s
#	      volGasFlow = kiln_data.inGasFlow / 1000.0, #m^3/s
#	      Tsolid = 420, #Tsolid K
#	      Tgas = 800, #Tgas K
#              Z0 = 0,
#              db=db)
#
#slice.T_wall = 369.2 #K
#slice.T_ext_shell = 0
#def f(x):
#    slice.T_ext_shell = x
#    return kiln.Q_w_ext(slice) - kiln.Q_w_sh(slice)
#
#slice.T_ext_shell = brenth(f, 298.15,  369.2, disp=True)
#
#print("T_ext=", slice.T_ext_shell)
#print("Q_w_ext=", kiln.Q_w_ext(slice))
#print("Q_w_sh=", kiln.Q_w_sh(slice))
#print("R_wall_shell=",kiln.R_wall_shell(slice))
#print("Q=",(slice.T_wall-slice.T_ext_shell)/kiln.R_wall_shell(slice))
#print("UA_shell_amb",1.0/(kiln.R_wall_shell(slice) + 1.0 / (1.0 / kiln.R_sh_ext_rd(slice) + 1.0 / kiln.R_sh_ext_cv(slice))))
#exit()

