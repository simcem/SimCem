import unittest

import simcem

class DatabaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = simcem.Database("free_database.xml")
        print(cls.db)
        
    def test_load(self):
        pass


###!/usr/bin/python3
##from pysimcem import *
##import unittest
##import time
##
##
##class TestSimCem(unittest.TestCase):
##    def setUp(self):
##        self.startTime = time.time()
##
##    def tearDown(self):
##        t = time.time() - self.startTime
##        print("%s: %.3fms" % (self.id(), t*1000))
##
##    def test_component_math(self):
##        A = Components({"CH4":1.0, "O2":1.0})
##        B = Components({"CO2":3.0, "O2":1.0, "NH3":0})
##        C = A+B
##        C.removeSmallComponents(0.01);
##        C = C * 10.0
##        C = 5 * C
##        C = C / 5.0
##        C_exact = Components({"CH4":10.0, "O2":20.0, "CO2":30.0})
##        self.assertEqual(C, C_exact)
##        self.assertEqual(C.N(), 60.0)
##        self.assertAlmostEqual(Components({"O2":1}).m(db), 0.032, places=5)
##
##    def test_phase_properties(self):
##        phase1 = ModelIdealGasTp(db, Components({"CO2":1.0}), 298.15, 1.01325e5) 
##        phase2 = ModelIdealGasTV(db, Components({"CO2":1.0}), 298.15, phase1.V())
##        self.assertAlmostEqual(phase1.Cp(),37.1354, places=2)
##        self.assertAlmostEqual(phase2.Cp(),37.1354, places=2)
##        self.assertAlmostEqual(phase1.H()/1000, -393510/1000, places=1)
##        self.assertAlmostEqual(phase2.H()/1000, -393510/1000, places=1)
##        self.assertAlmostEqual(phase1.S(), 213.7, places=1)
##        self.assertAlmostEqual(phase2.S(), 213.7, places=1)
##
##    def test_phase_properties_NASA_methane(self):
##        ##Test against the NASA online CEA tool (using same database)
##        phase1 = ModelIdealGasTp(db, Components({"CH4":1.0}), 350, 80e5)
##        phase2 = ModelIdealGasTV(db, Components({"CH4":1.0}), 350, phase1.V())
##        #CEA output is in kj/kg, need molar mass to convert
##        molar_mass = Components({"CH4":1.0}).M(db)
##        self.assertAlmostEqual(phase1.H()/molar_mass*1e-3, (-4531.31), places=1)
##        self.assertAlmostEqual(phase2.H()/molar_mass*1e-3, (-4531.31), places=1)
##        self.assertAlmostEqual(phase1.S()/molar_mass*1e-3, (9.7134), places=4)
##        self.assertAlmostEqual(phase2.S()/molar_mass*1e-3, (9.7134), places=4)
##        self.assertAlmostEqual(phase1.G()/molar_mass*1e-3, (-7931.02),places=1)
##        self.assertAlmostEqual(phase2.G()/molar_mass*1e-3, (-7931.02),places=1)
##        self.assertAlmostEqual(phase1.U()/molar_mass*1e-3, (-4712.71),places=1)
##        self.assertAlmostEqual(phase2.U()/molar_mass*1e-3, (-4712.71), places=1)
##        self.assertAlmostEqual(phase1.Cp()/molar_mass*1e-3, (2.3660), places=4)
##        self.assertAlmostEqual(phase2.Cp()/molar_mass*1e-3, (2.3660), places=4)
##        
##    def test_phase_properties_NASA_mixture(self):
##        #Try a mixture (equal weights)
##        n_CH4 = 1 / Components({"CH4":1.0}).M(db)
##        n_O2 = 1 / Components({"O2":1.0}).M(db)
##        phase1 = ModelIdealGasTp(db, Components({"CH4":n_CH4, "O2":n_O2}), 350, 80e5)
##        phase2 = ModelIdealGasTV(db, Components({"CH4":n_CH4, "O2":n_O2}), 350, phase1.V())
##        molar_mass = Components({"CH4":n_CH4, "O2":n_O2}).M(db)
##        
##        self.assertAlmostEqual(phase1.H()/molar_mass*1e-3, -2241.73, places=1)
##        self.assertAlmostEqual(phase2.H()/molar_mass*1e-3, -2241.73, places=1)
##        self.assertAlmostEqual(phase1.S()/molar_mass*1e-3, 7.8148, places=3)
##        self.assertAlmostEqual(phase2.S()/molar_mass*1e-3, 7.8148, places=3)
##        self.assertAlmostEqual(phase1.G()/molar_mass*1e-3, -4976.90, places=1)
##        self.assertAlmostEqual(phase2.G()/molar_mass*1e-3, -4976.90, places=1)
##        self.assertAlmostEqual(phase1.U()/molar_mass*1e-3, -2377.90, places=1)
##        self.assertAlmostEqual(phase2.U()/molar_mass*1e-3, -2377.90, places=1)
##        self.assertAlmostEqual(phase1.Cp()/molar_mass*1e-3, 1.6471, places=3)
##        self.assertAlmostEqual(phase2.Cp()/molar_mass*1e-3, 1.6471, places=3)
##
##    def test_class_overloading(self):
##        ###Test overloading of the classes
##        
##        class MyPhaseClassTV(ModelIdealGasTV):
##            def __init__(self, db, components, T, V):
##                self.volar_a_counter = 0
##                self.volar_u_counter = 0
##                self.CV_counter = 0
##                self.dpdV_counter = 0
##                self.dpdT_counter = 0
##                self.dpdNi_counter = 0
##                ModelIdealGasTV.__init__(self, db, components, T, V)
##                
##            def str(self):
##                #Don't use the "super" calling mechanism!
##                return "#" + ModelIdealGasTV.str(self)
##        
##            def volar_a(self, molID):
##                self.volar_a_counter = self.volar_a_counter + 1
##                return ModelIdealGasTV.volar_a(self, molID)
##        
##            def volar_u(self, molID):
##                self.volar_u_counter = self.volar_u_counter + 1
##                return ModelIdealGasTV.volar_u(self, molID)
##        
##            def CV(self):
##                self.CV_counter = self.CV_counter + 1
##                return ModelIdealGasTV.CV(self)
##        
##            def dpdV(self):
##                self.dpdV_counter = self.dpdV_counter + 1
##                return ModelIdealGasTV.dpdV(self)
##        
##            def dpdT(self):
##                self.dpdT_counter = self.dpdT_counter + 1
##                return ModelIdealGasTV.dpdT(self)
##        
##            def dpdNi(self, molID):
##                self.dpdNi_counter = self.dpdNi_counter + 1
##                return ModelIdealGasTV.dpdNi(self, molID)
##        
##        A = ModelIdealGasTV(db, Components({"N2":1.0}), 298.15, 10)
##        B = MyPhaseClassTV(db, Components({"N2":1.0}), 298.15, 10)
##        
##        self.assertEqual(A.str(),B.str()[1:])
##        self.assertEqual((A.volar_a("N2")),B.volar_a("N2"))
##        self.assertEqual(B.volar_a_counter,1)
##        self.assertEqual((A.chemPot("N2")), B.volar_a("N2"))
##        self.assertEqual(B.volar_a_counter, 2)
##        self.assertEqual((A.volar_u("N2")), B.volar_u("N2"))
##        self.assertEqual(B.volar_u_counter, 1)
##        self.assertEqual((A.dpdV()), B.dpdV())
##        self.assertEqual(B.dpdV_counter, 1)
##        self.assertEqual((A.dpdT()), B.dpdT())
##        self.assertEqual(B.dpdT_counter, 1)
##        self.assertEqual((A.dpdNi("N2")), B.dpdNi("N2"))
##        self.assertEqual(B.dpdNi_counter, 1)
##        self.assertEqual((A.CV()), B.CV())
##        self.assertEqual(B.CV_counter, 1)
##        self.assertEqual(B.dpdT_counter, 2)
##
##    #################################################
##    #################################################
##    #################################################
##    def test_Tp_water_flash(self):
##        # Check a simple Tp flash of water
##        # First check that its
##        for reactive in [True,False]:
##            sys = System(Objective_t.p, Objective_t.T, reactive)
##            A = ModelIdealGasTp(db, Components({"H2O":1}), 273.15+50, 1.01325e5)
##            B = ModelIncompressibleLiquid(db, Components({"H2O":1}), 273.15+50, 1.01325e5)
##            sys.append(A)
##            sys.append(B)
##            sys.equilibrate()
##            self.assertAlmostEqual(sys[0]["H2O"], 0, places=10)
##            self.assertAlmostEqual(sys[1]["H2O"], 2, places=10)
##            if not reactive:
##                self.assertAlmostEqual(-sys.lagrangians["H2O"]/sys[1].chemPot("H2O"), 1, places=9)
##
##    def test_pH_water_flash(self):
##        sys = System(Objective_t.p, Objective_t.H, reactive=False)
##        GasH = ModelIdealGasTp(db, Components({"H2O":1}), 273.15+100, 1.01325e5).H()
##        LiqH = ModelIncompressibleLiquid(db, Components({"H2O":1}), 273.15+100, 1.01325e5).H()
##
##        A = ModelIdealGasTp(db, Components({"H2O":0}), 273.15+100, 1.01325e5)
##        B = ModelIncompressibleLiquid(db, Components({"H2O":1}), 273.15+100, 1.01325e5)
##        sys.append(A)
##        sys.append(B)
##        sys.equilibrate(sys.p(), 0.5 * (GasH+LiqH))
##        self.assertAlmostEqual(sys[0]["H2O"], 0.5, places=2)
##        self.assertAlmostEqual(sys[1]["H2O"], 0.5, places=2)
##        lagrangian = sys[0].s("H2O")-sys[0].h("H2O")/sys[0].T()
##        self.assertAlmostEqual(sys.lagrangians["H2O"]/lagrangian, 1, places=7)
##            
##    def test_gas_Tp_combustion(self):
##        sys = System(Objective_t.p, Objective_t.T, True)
##        A = ModelIdealGasTp(db, Components({"N2":0.79, "O2":0.21, "CH4":0.105, "CO2":0.0, "CO":0.0, "H2O":0.0}), 2600, 1.01325e5)
##        sys.append(A)
##        #Test accessing the phase via the PhaseBase interface.
##        # System is also a MultiPhase container, and any access to its
##        # contained elements provides a shared_ptr<PhaseBase> type.
##        a_base = sys[0]
##        self.assertEqual(a_base.T(), 2600)
##        self.assertEqual(a_base.p(), 1.01325e5)
##        sys.equilibrate()
##        testResult = ModelIdealGasTp(db, Components({"CH4":8.688777086942181e-12, "CO":0.03453317988062815, "CO2":0.07046682011068307, "H2O":0.2099999999826224, "N2":0.79, "O2":0.01726658995769161}), 2600, 1.01325e5);
##        result = testResult - A
##        result.removeSmallComponents(1e-7)
##        self.assertEqual(len(result), 0, msg="remainder "+str(result))
##
##    def test_gas_TV_combustion(self):
##        sys = System(Objective_t.p, Objective_t.T, True)
##        V = (0.79+0.21+0.105) * 8.3145 * 2600 / 1.01325e5
##        A = ModelIdealGasTV(db, Components({"N2":0.79, "O2":0.21, "CH4":0.105, "CO2":0.0, "CO":0.0, "H2O":0.0}), 2600, V)
##        sys.append(A)
##        sys.equilibrate()
##        testResult = ModelIdealGasTp(db, Components({"CH4":8.678512425246362e-12, "CO":0.03453322192945054, "CO2":0.07046677806187095, "H2O":0.209999999982643, "N2":0.79, "O2":0.01726661098208228}), 2600, 1.01325e5);
##        result = testResult - A
##        result.removeSmallComponents(1e-7)
##        self.assertEqual(len(result), 0, msg="remainder "+str(result))
##                
##    def test_solidgas_Tp_combustion(self):
##        #Solid combustion
##        sys = System(Objective_t.p, Objective_t.T, True)
##        A = ModelIdealGasTp(db, Components({"O2":1.0, "CO2":0, "CO":0}), 2600, 1.01325e5)
##        B = ModelIncompressible(db, Components({"C:graphite":1.0, "CaO:Crystal":0}), 2600, 1.01325e5, "solid", True)
##        sys.append(A)
##        sys.append(B)
##        sys.equilibrate()
##        testResult = ModelIdealGasTp(db, Components({"CO":0.1760714215970302, "CO2":0.823928578398245, "O2":0.08803571080323989}), 2600, 1.01325e5);
##        result = testResult - A
##        result.removeSmallComponents(1e-6)
##        self.assertEqual(len(result), 0, msg="remainder "+str(result))
##
##    def test_solidgas_TV_combustion(self):
##        sys = System(Objective_t.p, Objective_t.T, True)
##        V = 1.0 * db.R * 2600 / 1.01325e5
##        A = ModelIdealGasTV(db, Components({"O2":1.0, "CO2":0, "CO":0}), 2600, V)
##        B = ModelIncompressible(db, Components({"C:graphite":1.0, "CaO:Crystal":0}), 2600, 1.01325e5, "solid", True)
##        sys.append(A)
##        sys.append(B)
##        sys.equilibrate()
##        testResult = ModelIdealGasTp(db, Components({"CO":0.1760715474833369, "CO2":0.823928452512056, "O2":0.08803577374627462}), 2600, 1.01325e5);
##        result = testResult - A
##        result.removeSmallComponents(1e-6)
##        self.assertEqual(len(result), 0, msg="remainder "+str(result))
##
##    def test_gas_adiabatic_Tp_combustion(self):
##        #Adiabatic combustion
##        sys = System(Objective_t.p, Objective_t.H, True)
##        A = ModelIdealGasTp(db, Components({"N2":0.79, "O2":0.21, "CH4":0.105, "CO2":0.0, "CO":0.0, "H2O":0.0}), 300, 1.01325e5)
##        sys.append(A)
##        sys.equilibrate()
##        testResult = ModelIdealGasTp(db, Components({"CH4":6.294094972974197e-12, "CO":0.01147139046161158, "CO2":0.09352860953209433, "H2O":0.2099999999874118, "N2":0.79, "O2":0.005735695243393977}), 2258.08, 1.01325e5);
##        result = testResult - A
##        result.removeSmallComponents(1e-7)
##        self.assertEqual(len(result), 0, msg="remainder "+str(result))
##        self.assertAlmostEqual(A.T(), testResult.T(), places=1)
##
##    def test_gas_adiabatic_TV_combustion(self):
##        sys = System(Objective_t.p, Objective_t.H, True)
##        V = (0.79+0.21+0.105) * 8.3145 * 300 / 1.01325e5
##        sys.append(ModelIdealGasTV(db, Components({"N2":0.79, "O2":0.21, "CH4":0.105, "CO2":0.0, "CO":0.0, "H2O":0.0}), 300, V))
##        sys.equilibrate()
##        testResult = ModelIdealGasTp(db, Components({"CH4":6.294094972974197e-12, "CO":0.01147139046161158, "CO2":0.09352860953209433, "H2O":0.2099999999874118, "N2":0.79, "O2":0.005735695243393977}), 2258.08, 1.01325e5);
##        result = testResult - sys[0]
##        result.removeSmallComponents(2e-6)
##        self.assertEqual(len(result), 0, msg="remainder "+str(result))
##        self.assertAlmostEqual(sys[0].T(), testResult.T(), places=1)
##
##    def test_VT_Tp_combustion(self):
##        #Isothermal-Isochoric combustion
##        sys = System(Objective_t.V, Objective_t.T, True)
##        A = ModelIdealGasTp(db, Components({"N2":0.79, "O2":0.21, "H2":0.42, "H2O":0.0}), 2000, 1.01325e5)
##        origV = A.V()
##        orig_elements = A.components.elements(db)
##        sys.append(A)
##        sys.equilibrate()
##        testResult = ModelIdealGasTp(db, Components({"H2":0.003476625431155716, "H2O":0.4165233745688441, "N2":0.79, "O2":0.001738312715577942}), 2000, 0.864644e5)
##        result = testResult - A
##        result.removeSmallComponents(1e-5)
##        self.assertEqual(len(result), 0, msg="remainder "+str(result))
##        self.assertAlmostEqual(A.V(), origV, places=4)
##
##    def test_getProperty(self):
##        self.assertEqual(db.getComponent("C6H12,cyclo-").getProperties()[Property_t.pc][0].value, 4079999.9999999995)
##
    
if __name__ == '__main__':
    unittest.main()
