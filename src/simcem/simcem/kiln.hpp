#pragma once
#include <simcem/simcem.hpp>
#include <simcem/sundials.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/constants/constants.hpp>
#include <memory>
#include <numeric>

/*
  Notes:
  Theo's guide

  There are three PyKiln folders:  
  - PyKiln - Output from each test
  - PyKiln2 - Output from each test
  - Optimize_Heat_Model - Theo's work on optimising the heat transfer coefficients

  You need to install assimulo to use it. Do this by downloading the
  source code, replacing all StrictVersion by LooseVersion in
  setup.py, then build/install.

  You also need to install PyChemEng. The folders in dropbox are not
  installable (old versions) It seems that PyChemEngLocal is the
  newest version, followed by PyChemEng-Marcus, then PyChemEng. You
  need to install from git@dynamomd.org:PyChemEng.

  You also need pyopt!

  Each PyKiln folder has the following structure:
  - Kiln_Classes.py - All pilot kiln trial data.
  - KilnSolverXXXX.py - The actual workhorse for fitting data, all copies of KilnSolverXXX.py
  - Kiln_Plotter.py - Plotting of data files from KilnSolver
  
  run PyKiln2/KilnSolverAssimulo.py, then run Kiln_Plotter.py
*/
namespace simcem {
  /*! \brief Heat transfer correlations and expressions. */
  namespace HT {
    /*! \brief Conductive heat transfer. */
    namespace conduction {
      /*! \brief Conduction resistance for a cylindrical shell of unit length.*/
      double R_cyl(const double Rout, const double Rin, const double k) {
	return std::log(Rout / Rin) / (2 * Database::pi * k);
      }

      /*! \brief The Maxwell model for combining thermal conductivities. 
	
	This requires that one phase is identified as "continuous" and
	another phase is "discrete". Then the volume fraction of
	discontinuous phase is requred.
      */
      double maxwell(const double kcont, const double kdisc, const double discfrac) {
	return kcont * (2*kcont + kdisc + 2 * discfrac * (kdisc - kcont)) / (2 * kcont + kdisc + discfrac * (kcont - kdisc));
      }
    }
    
    /*! \brief Natural convection heat transfer. */
    namespace natConv {
      /*! \brief The overall convective heat transfer from smooth
	circular cylinders.
	
	Taken from Heat Transfer, J.P Holman pg 334, VT Morgan.

	A big collection of available expressions is in Boetcher_2014. 
      */
      double Nu_horiz_cyl(const double Gr, const double Pr) {
	const double Ra = Gr * Pr;

	//A wide range expression (Eq.7-36 in Holman, Eq 5-35 in Perrys)
	//	if ((1e-5 <= Ra) && (Ra <= 1e12))
	//	  return std::pow(0.60 + 0.387*std::pow(Ra/std::pow(1+std::pow(0.559/Pr, 9.0/16), 16/9.0), 1.0/6), 2);

	//The table from Holman favouring the "preferred" expressions where possible
	if (Ra < 1e-10)
	  return 0.4;
	if (Ra < 1e-2)
	  return 0.675 * std::pow(Ra, 0.058);
	if (Ra < 1e2)
	  return 1.02 * std::pow(Ra, 0.148);
	if (Ra < 1e4)
	  return 0.850 * std::pow(Ra, 0.188);
	if (Ra < 1e7)
	  return 0.480 * std::pow(Ra, 1.0/4);
	if (Ra < 1e12)
	  return 0.125 * std::cbrt(Ra);
       	
	stator_throw() << "Ra=" << Ra << " over upper range.";
      }
    }

    /*!\brief Radiative heat transfer. */
    namespace rad {

      /*! \brief Simple radiative heat transfer coefficient.
	
	The Steffan-Boltzmann coefficient 
      */
      double h_rad_simple(const double T1, const double T2, const double emissivity,
			  const double steffanBoltzmannConstant = 5.670367e-8) {
	return emissivity *
	  steffanBoltzmannConstant
	  * (T1 * T1 + T2 * T2) * (T1 + T2);
      }

      
      /*! \brief Hottel calculation for the emissivity of gases
	containing CO2 and H2O.
	
	This is taken from Perry's Chemical Engineering Handbook,
	chapter 5. This assumes that the primary radiative
	properties of the gas are derived from the H2O and CO2
	components only. It also uses bilinear
	interpolation/extrapolation.

	Errors up to 20% are expected. The temperature argument is
	only provided for absorptivity calculations (which should be
	carried out via the Hottel_absorptivity function), and if it
	is left at the default of 0 it is taken from the gas phase
	automatically.
      */
      double Hottel_emissivity(const simcem::Model& gas, const double pathLength, double Tgas = 0) {
	if (Tgas == 0) Tgas = gas.T();
	//Fraction of H2O versus CO2+H2O
	const double pw_frac = gas["H2O"] / (gas["H2O"]+gas["CO2"]);

	//Partial pressure of H2O and CO2 in atms
	const double p = (gas["H2O"]+gas["CO2"]) / gas.N() * gas.p() / 1.01325e5;
	const double pL = p * pathLength;
	if (pL > 10)
	  stator_throw() << "pL out of range " << pL;

	if (pL < 0.005)
	  return 0; //This is an approximation!
	  //stator_throw() << "pL out of range " << pL;
	
	const double logpL = std::log10(pL);

	//Perform bilinear interpolation where x is T, and y is pw_frac.
	const double x_range[] = {1000, 1500, 2000};
	const double y_range[] = {0, 1.0/3.0, 1.0/2.0, 2.0/3.0, 3.0/4.0, 1.0};

	//Calculation of the indices of the data points are around the sample point
	size_t x_idx = 0;
	if (Tgas > x_range[1]) x_idx = 1;
	
	size_t y_idx = 0;
	if (pw_frac > y_range[1]) y_idx = 1;
	if (pw_frac > y_range[2]) y_idx = 2;
	if (pw_frac > y_range[3]) y_idx = 3;
	if (pw_frac > y_range[4]) y_idx = 4;
	
	static const double coeffs[] = { 2.2661, 0.1742, -0.0390,  0.0040,
					 2.3954, 0.2203, -0.0433,  0.00562,
					 2.4104, 0.2602, -0.0651, -0.00155,
					 2.5754, 0.2792, -0.0648,  0.0017,
					 2.6451, 0.3418, -0.0685, -0.0043,
					 2.6504, 0.4279, -0.0674, -0.0120,
					 2.6090, 0.2799, -0.0745, -0.0006,
					 2.6862, 0.3450, -0.0816, -0.0039,
					 2.7029, 0.4440, -0.0859, -0.0135,
					 2.6367, 0.2723, -0.0804,  0.0030,
					 2.7178, 0.3386, -0.0990, -0.0030,
					 2.7482, 0.4464, -0.1086, -0.0139,
					 2.6432, 0.2715, -0.0816,  0.0052,
					 2.7257, 0.3355, -0.0981,  0.0045,
					 2.7592, 0.4372, -0.1122, -0.0065,
					 2.5995, 0.3015, -0.0961,  0.0119,
					 2.7083, 0.3969, -0.1309,  0.00123,
					 2.7709, 0.5099, -0.1646, -0.0165};

	//A routine to calculate the offset for the coefficients
	auto C = [&](size_t x_idx, size_t y_idx) { return coeffs + 4 * (x_idx + 3 * y_idx); };

	//A calculation of the function at the index point
	auto f = [&](const double* C) { return C[0] + logpL * (C[1] + logpL * (C[2] + logpL * C[3])); };
	
	const double x = (Tgas - x_range[x_idx]) / (x_range[x_idx+1] - x_range[x_idx]);
	const double y = (pw_frac   - y_range[y_idx]) / (y_range[y_idx+1] - y_range[y_idx]);

	double log10eT = 
	  f(C(x_idx, y_idx)) * (1 - x) * (1 - y)
	  + f(C(x_idx + 1, y_idx)) * x * (1 - y)
	  + f(C(x_idx, y_idx + 1)) * (1 - x) * y
	  + f(C(x_idx + 1, y_idx + 1)) * x * y;
	
	return std::pow(10, log10eT) / Tgas;
      }

      /* \brief Hottel calculation for the absorptivity of gases
	 containing CO2 and H2O.

	 This assumes ideal gas behaviour!
      */
      
      double Hottel_absorptivity(simcem::Model& gas, double pathLength, double Tsurface) {
	//Adjust pathlength, such that pL is scaled (this is a
	//horrible hack to adjust the pressure in the emissivity
	//calculations according to the ideal gas law!)
	const double L = pathLength * Tsurface / gas.T();
	const double e = Hottel_emissivity(gas, L, Tsurface);
	return e * std::sqrt(gas.T() / Tsurface);
      }

      /*! \brief Two-body radiation problem in third region.
	
	\image html twobodygas_radiation_network.svg
	
	We assume all resistances and the voltages \f$V_a\f$, \f$V_b\f$, and \f$V_c\f$ are known.
	The total currents at the nodes \f$\alpha\f$ and \f$\beta\f$ must sum to zero, thus:
	
	\f[
	(V_a - V_\alpha)R^{-1}_1 + (V_\beta - V_\alpha)R^{-1}_2 + (V_c - V_\alpha)R^{-1}_3 = 0
	\f]
	\f[
	(V_b - V_\beta)R^{-1}_4 + (V_\alpha - V_\beta)R^{-1}_2 + (V_c - V_\beta)R^{-1}_5 = 0
	\f]
	
	Solving for \f$V_\beta\f$:
	
	\f[
	V_\beta = \frac{m_2\,c_2+c_1}{1-m_1\,m_2}
	\f]
	\f[
	V_\alpha = m_2\, V_\beta + c_2
	\f]
	
	where:
	\f[
	m_1 = \frac{1}{R_2^{-1}}\left(R_1^{-1} + R_2^{-1} + R_3^{-1}\right)
	\f]
	\f[
	c_1 = -V_a\frac{R_1^{-1}}{R_2^{-1}} -V_c\frac{R_3^{-1}}{R_2^{-1}}
	\f]
	
	\f[
	m_2 = \frac{1}{R_2^{-1}}\left(R_4^{-1} + R_2^{-1} + R_5^{-1}\right)
	\f]
	\f[
	c_2 = -V_b\frac{R_4^{-1}}{R_2^{-1}} -V_c\frac{R_5^{-1}}{R_2^{-1}}
	\f]

	The currents are returned as follows:
	\f[
	I_a = R_1^{-1}(V_\alpha-V_a) \qquad I_b = R_4^{-1}(V_\beta-V_b) \qquad I_c = R_3^{-1}(V_\alpha-V_c) + R_5^{-1}(V_\beta-V_c)
	\f]
	
	As a test, we should have \f$I_a+I_b+I_c=0\f$.
      */
      std::array<double, 3> twobodygas_network(const std::array<double, 3> V, const std::array<double, 5> Rinv) {
	const double m1 = (Rinv[0]+Rinv[1]+Rinv[2]) / Rinv[1];
	const double m2 = (Rinv[3]+Rinv[1]+Rinv[4]) / Rinv[1];
	const double c1 = -V[0]*Rinv[0]/Rinv[1] - V[2] * Rinv[2] / Rinv[1];
	const double c2 = -V[1]*Rinv[3]/Rinv[1] - V[2] * Rinv[4] / Rinv[1];
	const double Vbeta = (m1 * c2 + c1) / (1 - m1 * m2);
	const double Valpha = m2 * Vbeta + c2;

	return std::array<double, 3>{
	  {Rinv[0] * (Valpha - V[0]),
	      Rinv[3] * (Vbeta - V[1]),
	      Rinv[2] * (Valpha - V[2]) + Rinv[4] * (Vbeta - V[2])
	      }
	};
      }
    }

    /*! \brief Specialised expressions for heat transfer in kilns. */
    namespace kiln {
      /*! \brief Nusselt number correlation for heat transfer between
	a cylindrical wall and a granular bed of material which
	covers it.

	This is taken from Li_etal_2005. The Nusselt number uses the
	gas thermal conductivity and particle diameter as
	normalisation factors.

	\param Peclet The Peclet number for the bed (see Pe_bed_wall_Li_etal_2005).
	\param chi A fitting parameter, which appears to be around 0.1, but is in the range 0.096-0.198.
      */
      double Nu_bed_wall_Li_etal_2005(const double Peclet, const double chi = 0.1) {
	return 1 / (chi + 0.5 * std::sqrt(Database::pi / Peclet));
      }

      
      /*! \brief Peclet number for Heat transfer between a cylindrical
	wall and a granular bed of material which covers it.

	This is taken from Li_etal_2005.

	\param dp Particle diameter \f$d_p\f$.
	\param kg Gas thermal conductivity \f$k_g\f$.
	\param kb Solid bed thermal conductivity \f$k_b\f$.
	\param vol_Cp Solid bed volumetric heat capacity, \f$\rho_b\,c_{p,b}\f$.
	\param centralangle Central angle of bed, \f$\theta\f$.
	\param angularvel Angular velocity of cylinder, \f$\omega\f$.
      */
      double Pe_bed_wall_Li_etal_2005(const double dp, const double kg, const double kb,
				      const double vol_Cp, const double centralangle, const double angularvel ) {
	return std::pow(dp/kg, 2) * vol_Cp * kb * angularvel / centralangle;
      }

      /*! \brief Heat transfer coefficient between a wall and a
	granular bed of material which covers it.

	This is taken from Li_etal_2005.

	\param dp Particle diameter \f$d_p\f$.
	\param kg Gas thermal conductivity \f$k_g\f$.
	\param kb Solid bed thermal conductivity \f$k_b\f$.
	\param vol_Cp Solid bed volumetric heat capacity, \f$\rho_b\,c_{p,b}\f$.
	\param centralangle Central angle of bed, \f$\theta\f$.
	\param angularvel Angular velocity of cylinder, \f$\omega\f$.
	\param chi A fitting parameter, which appears to be around 0.1, but is in the range 0.096-0.198.
      */
      double h_bed_wall_Li_etal_2005(const double dp, const double kg, const double kb,
				     const double vol_Cp, const double centralangle, const double angularvel,
				     const double chi = 0.1) {
	const double Pe = Pe_bed_wall_Li_etal_2005(dp, kg, kb, vol_Cp, centralangle, angularvel);
	const double Nu = Nu_bed_wall_Li_etal_2005(Pe, chi);
	return Nu * kg / dp;
      }

      /*! \brief Gas-wall Nusselt number for rotary kilns.
	
	\f[
	\text{Nu}_{g-ew} = 1.54\,\text{Re}_g^{0.575}\,\text{Re}_\omega^{-0.292}
	\f]
	
	This is taken from Eq.9 of Li_etal_2005. This requires
	\f$1600<\text{Re}_g<7800\f$ and \f$20<\text{Re}_\omega<800\f$.
	
	\param Re_g Reynolds number for the gas, \f$\text{Re}_g\f$, using the normal hydraulic diameter.
	\param Re_w Angular rotation Reynolds number, \f$\text{Re}_\omega\f$.
      */
      double Nu_gas_wall_Li_etal_2005(const double Re_g, const double Re_w) {
#ifndef NDEBUG
	//if (!((1600 < Re_g) && (Re_g < 7800)))
	//stator_throw() << "Out of range for Re_g " << Re_g;

	//if (!((20 < Re_w) && (Re_w < 800)))
	//stator_throw() << "Out of range for Re_w " << Re_w;
#endif

	return 1.54 * std::pow(Re_g, 0.575) * std::pow(Re_w, -0.292);
      }

      /*! \brief Gas-bed Nusselt number for rotary kilns.

	\f[
	\text{Nu}_{g-eb} = 0.46\,\text{Re}_g^{0.535}\,\text{Re}_\omega^{0.104}\,\eta^{-0.341}
	\f]

	This is taken from Eq.8 of Li_etal_2005. This requires
	\f$1600<\text{Re}_g<7800\f$ and \f$20<\text{Re}_\omega<800\f$.
	
	\param Re_g Reynolds number for the gas, \f$\text{Re}_g\f$, using the normal hydraulic diameter.
	\param Re_w Angular rotation Reynolds number, \f$\text{Re}_\omega\f$.
	\param bed_frac Fraction of kiln filled with the bed, \f$\eta\f$.
      */
      double Nu_gas_bed_Li_etal_2005(const double Re_g, const double Re_w, const double bed_frac) {
#ifndef NDEBUG
//	if (!((1600 < Re_g) && (Re_g < 7800)))
//	  stator_throw() << "Out of range for Re_g " << Re_g;
//
//	if (!((20 < Re_w) && (Re_w < 800)))
//	  stator_throw() << "Out of range for Re_w " << Re_w;
#endif
	return 0.46 * std::pow(Re_g, 0.535) * std::pow(Re_w, 0.104) * std::pow(bed_frac, -0.341);
      }

      
      /*! \brief Correlation for mean beam length in a rotary kiln.
	
	This is taken from Eq.50 of Gorog_etal_1981.

	\param R Kiln radius.
	\param h Bed height (at its deepest/central point);
      */
      double pathLength(const double R, const double h) {
	return 2 * R * 0.95 * (1 - h / (2 * R));
      }
    }
  }

  
  
  namespace kiln {
    struct Slice {
      Slice(simcem::shared_ptr<simcem::ModelIdealGasTp> gas,
	    simcem::shared_ptr<simcem::ModelIncompressible> solid,
	    double Z = 0
	    ):
	_Z(Z)
      {
	//Here we copy the phases to ensure the slice has a unique copy
	_gas = simcem::shared_ptr<simcem::ModelIdealGasTp>(new simcem::ModelIdealGasTp(*gas));
	_solid = simcem::shared_ptr<simcem::ModelIncompressible>(new simcem::ModelIncompressible(*solid));

	_T_ext_shell = 0.75*298.15 + 0.25 * _gas->T();
	_T_wall = 0.25*298.15 + 0.75 * _gas->T();
      }

      Slice(const Slice& s):
	_gas_k(s._gas_k),
	_gas_visc(s._gas_visc),
	_gas_density(s._gas_density),
	_gas_vel(s._gas_vel),
	_solid_k(s._solid_k),
	_solid_Cp(s._solid_Cp),
	_T_ext_shell(s._T_ext_shell),
	_T_wall(s._T_wall),
	_Z(s._Z)
      {
	_gas = simcem::shared_ptr<simcem::ModelIdealGasTp>(new simcem::ModelIdealGasTp(*s._gas));
	_solid = simcem::shared_ptr<simcem::ModelIncompressible>(new simcem::ModelIncompressible(*s._solid));
      }

      Slice(const simcem::Components solid,
	    const double volAirFlow,
	    const double volGasFlow,
	    const double volO2Flow,
	    const double volSO2Flow,
	    const double Tsolid,
	    const double Tgas,
	    const double Z0,
	    const simcem::shared_ptr<simcem::Database> db):
	_Z(Z0)
      {
	const simcem::ModelIdealGasTp air(db, Database::getDryAir(), 298.15, 1.01325e5);
	const simcem::Components inletair = Database::getDryAir() * ((volAirFlow) / air.V());
    
	const simcem::ModelIdealGasTp oxy(db, Components({{"O2",1}}), 298.15, 1.01325e5);
	const simcem::Components inletoxy = Components({{"O2",1}}) * ((volO2Flow) / oxy.V());
	
	const simcem::ModelIdealGasTp so2(db, Components({{"SO2",1}}), 298.15, 1.01325e5);
	const simcem::Components inletso2 = Components({{"SO2",1}}) * ((volSO2Flow) / oxy.V());

	const simcem::ModelIdealGasTp fuel(db, {{"CH4",100.0}}, 298.15, 1.01325e5);
	const simcem::Components inletfuel = Components(fuel) * (volGasFlow / fuel.V());

	const simcem::Components combustion_outputs({{"H2O",0}, {"CO2", 0}});

	_gas = simcem::shared_ptr<simcem::ModelIdealGasTp>(new simcem::ModelIdealGasTp(db, inletfuel + inletair + inletoxy + inletso2 + combustion_outputs, 298.15, 1.01325e5));
	simcem::System sys(simcem::Objective_t::p, simcem::Objective_t::H, true);
	sys.push_back(_gas);

	try {
	  sys.equilibrate();
	} catch (std::exception& e) {
	  std::cout << "Failed to equilibrate initial gas " << _gas->str() << std::endl;
	}

	_solid = simcem::shared_ptr<simcem::ModelIncompressible>(new simcem::ModelIncompressible(db, solid, 298.15, 1.01325e5, "solid"));

	//Value of Tgas=0 means retain the adiabatic temperature
	if (Tgas != 0)
	  _gas->set(simcem::Objective_t::T, Tgas, simcem::Objective_t::p);
	
	_solid->set(simcem::Objective_t::T, Tsolid, simcem::Objective_t::p);
    
	_T_ext_shell = 0.75*298.15 + 0.25 * _gas->T();
	_T_wall = 0.25*298.15 + 0.75 * _gas->T();
      }
      
      simcem::shared_ptr<simcem::ModelIdealGasTp> _gas;
      simcem::shared_ptr<simcem::ModelIncompressible> _solid;

      /*! \brief Precompute a few frequently used values for the slice. */
      void compute_properties(const double gas_area, const double solid_k) {
	//Gas properties
	std::tie(_gas_visc, _gas_k) = simcem::trans::NASA_transport(_gas);
	_gas_density = _gas->M() / _gas->V();

	//Solid properties
	_solid_Cp = _solid->Cp() / _solid->M();
	_gas_vel = _gas->V() / gas_area;
	_solid_k = solid_k;
      }
	
      //Cached values of the computed properties of the phases
      double _gas_k;
      double _gas_visc;
      double _gas_density;
      double _gas_vel;
      double _solid_k;
      double _solid_Cp;
      double _T_ext_shell;
      double _T_wall;
      double _Z;
    };
        
    /*! \brief Abstract base class for models describing the solid bed dynamics.*/
    class BedModel {
    public:
      virtual double bedHeight(const Slice& slice) const = 0;
    };

    class ConstantHeightBed : public BedModel {
    public:
      ConstantHeightBed(double height):
	_height(height)
      {}

      virtual double bedHeight(const Slice& slice) const { return _height; }

      static simcem::shared_ptr<BedModel> fromFillingFraction(double frac, double kiln_radius) {
	const double pi = boost::math::constants::pi<double>();
	
	std::pair<double, double> central_angle_limits
	  = boost::math::tools::bisect([pi, frac](double x) {return x - std::sin(x) - 2 * pi * frac;},
				       0.0, 2 * pi,
				       boost::math::tools::eps_tolerance<double>(10));

	double central_angle = (central_angle_limits.first + central_angle_limits.second) / 2;

	return simcem::shared_ptr<BedModel>
	  (new ConstantHeightBed(kiln_radius * (1.0 - std::cos(central_angle/2))));
      }
      
    protected:
      double _height;
    };


    /*! \brief A process model for a rotary kiln.
     
      \image html kiln_balance.svg
      
      Performing a differential balance of mass over a phase \f$\alpha\f$:
      \f[
      \Delta z\,A_\alpha\left(\rho_\alpha(z,\,t+\Delta t)-\rho_\alpha(z,\,t)\right)=  \Delta t\,A_\alpha\left(\dot{m}_\alpha\left(z\right)-\dot{m}_\alpha\left(z+\Delta z\right)\right) + \Delta t\,\Delta z\,\dot{m}_{\beta\to\alpha}
      \f]

      where \f$\dot{m}_{\beta\to\alpha}\f$ is the rate of mass
      transfer between the two phases per unit length of kiln and
      conservation of mass requires
      \f$\dot{m}_{\beta\to\alpha}=-\dot{m}_{\alpha\to\beta}\f$. Dividing
      by \f$\Delta t\,\Delta z\f$ and taking the limit as they go to
      zero:

      \f[
      \frac{{\rm d}\,A_\alpha\,\rho_\alpha}{{\rm d}t}= -\frac{{\rm d} A_\alpha\,\dot{m}_\alpha}{{\rm d}z} + \dot{m}_{\beta\to\alpha}
      \f]

      For enthalpy:
      \f[
      \Delta z\,A_\alpha\left(\rho_\alpha(z,\,t+\Delta t)\,h_\alpha(z,\,t+\Delta t)-\rho_\alpha(z,\,t)\,h_\alpha(z,\,t)\right)=  \Delta t\,A_\alpha\left(\dot{m}_\alpha\left(z\right)\,h_\alpha(z,\,t)-\dot{m}_\alpha\left(z+\Delta z\right)\,h_\alpha(z+\Delta z,\,t)\right) + \Delta t\,\Delta z\,Q_{\to\alpha}
      \f]
      where \f$Q_{\to\alpha}\f$ is the heat transfer rate to phase \f$\alpha\f$ per unit length of kiln. Again, a differential form results:
      \f[
      \frac{{\rm d}\,A_\alpha\,\rho_\alpha\,h_\alpha}{{\rm d}t}= -\frac{{\rm d}\,A_\alpha\,\dot{m}_\alpha\,h_\alpha}{{\rm d}z} + Q_{\to\alpha}

      
      \f]
    */
    class Kiln {
    public:
      struct Layer {
	std::string _material;
	double _thickness;
	sym::Expr _k;
      };
    
      Kiln(double RPM, double innerRadius, double length,
	   double particle_diam,
	   double shell_emissivity,
	   double bed_emissivity,
	   double wall_emissivity,
	   double solid_density,
	   double bed_void_frac,
	   sym::Expr solid_k,
	   simcem::shared_ptr<Database> db
	   ):
	_RPM(RPM), _innerRadius(innerRadius), _length(length),
	_particle_diam(particle_diam),
	_shell_emissivity(shell_emissivity),
	_bed_emissivity(bed_emissivity),
	_wall_emissivity(wall_emissivity),
	_solid_density(solid_density),
	_bed_void_frac(bed_void_frac),
	_solid_k(solid_k),
	_ambient(new simcem::ModelIdealGasTp(db, Database::getDryAir(), 298.15, 1.01325e5))
      {}

      void setFixedHeightBedModel(const double solidLoading) {
	_bedModel = simcem::kiln::ConstantHeightBed::fromFillingFraction(solidLoading, _innerRadius);
      }
      
      void add_layer(std::string material, double thickness, sym::Expr k) {
	_layers.push_back(Layer{material, thickness, k});
      }
      
      double solid_k(const double T) {
	return sym::fast_sub(_solid_k, sym::Var<sym::vidx<'T'>>() = T);
      }

      double innerRadius() const { return _innerRadius; }
    
      double outerRadius() const {
	return _innerRadius + std::accumulate(_layers.begin(), _layers.end(), 0.0, [](double a, Layer b){ return a + b._thickness; });
      }

      /*! \brief Conductive heat transfer resistance in the kiln wall. */
      double R_wall_shell(const Slice& slice) const {
	double Ri = _innerRadius;
	double Rtotal = 0;
	for (const auto& layer : _layers) {
	  const double Ro = Ri + layer._thickness;
	  const double k = sym::fast_sub(layer._k, sym::Var<sym::vidx<'T'>>() = slice._T_wall);
	  Rtotal += std::log(Ro / Ri) / (2 * simcem::Database::pi * k);
	  Ri = Ro;
	};

	return Rtotal;
      }
      
      const BedModel& bedModel() const { return *_bedModel; }

      const std::vector<Slice>& slices() const { return _slices; }
      std::vector<Slice>& slices() { return _slices; }
      
      /*! \brief The chord length of the bed.
	
	The chord length can be calculated in two ways:
	\f[
	L_{chord}= 2\,R\,\sin(\theta/2) = 2\sqrt{h(2*R-h)}
	\f]

	The latter is used here.
      */
      double chordLength(const Slice& slice) const {
	const double h = _bedModel->bedHeight(slice);
	return 2 * std::sqrt(h * (2 * _innerRadius - h));
      }

      /*! \brief The central angle is the radians of the cylinder
	covered by the bed.
      */
      double centralAngle(const Slice& slice) const {
	return 2 * std::asin(chordLength(slice) / 2 / _innerRadius);
      }

      double bedArea(const Slice& slice) const {
	const double cAngle = centralAngle(slice);
	return 0.5 * _innerRadius * _innerRadius * (cAngle - std::sin(cAngle));
      }

      double kilnArea() const {
	return simcem::Database::pi * _innerRadius * _innerRadius;
      }

      double gasArea(const Slice& slice) const {
	return kilnArea() - bedArea(slice);
      }

      double gasPerimeter(const Slice& slice) const {
	return (2 * simcem::Database::pi - centralAngle(slice)) * _innerRadius + chordLength(slice);
      }
      
      double hydraulicDiameter(const Slice& slice) const {
	return 4 * gasArea(slice) / gasPerimeter(slice);
      }

      double exposedWallPerimeter(const Slice& slice) const {
	const double pi = boost::math::constants::pi<double>();
	const double cAngle = centralAngle(slice);
	return (2 * pi - cAngle) * _innerRadius;
      }

      double coveredWallPerimeter(const Slice& slice) const {
	const double cAngle = centralAngle(slice);
	return cAngle * _innerRadius;
      }
      
      double angular_vel() const {
	return _RPM / 60.0 * 2 * simcem::Database::pi;
      }
      
      /*! \brief Heat transfer coefficient between the bed and the wall it covers, \f$h_{cw-s}\f$.
       */
      double h_bed_wall(const Slice& slice) const {
	const double Cp_per_vol = slice._solid_Cp * _solid_density;
	const double kb = HT::conduction::maxwell(slice._gas_k, slice._solid_k, _bed_void_frac);
	return simcem::HT::kiln::h_bed_wall_Li_etal_2005(_particle_diam, slice._gas_k, kb,
							 Cp_per_vol, centralAngle(slice), angular_vel());
      }

      /*! \brief \f$Q_{w\to s}^{cd}\f$ */
      double Q_wall_bed_cd(const Slice& slice) const {
	return h_bed_wall(slice) * coveredWallPerimeter(slice) * (slice._T_wall - slice._solid->T());
      }
      
      /*! \brief Free gas Reynolds number, \f$\text{Re}_g\f$.
	
	This is the Reynolds number for the gas above the solid bed.
	
	\f[
	\text{Re}_g = \frac{\rho_g\,v_g\,D_e}{\mu_g}
	\f]

	where \f$D_e\f$ is the hydraulic diameter.
      */
      double Re_g(const Slice& slice) const {
	const double De = hydraulicDiameter(slice);
	return slice._gas_density * slice._gas_vel * De / slice._gas_visc;
      }
      
      /*! \brief Rotational gas Reynolds number, \f$\text{Re}_\omega\f$.
	
	This is the Reynolds number for the gas within the bed itself.
	
	\f[
	\text{Re}_\omega = \frac{\rho_g\,\omega\,D_e^2}{\mu_g}
	\f]

	where \f$D_e\f$ is the hydraulic diameter. An example of this
	definition is given just below Eq.9 in Li_etal_2005.
      */
      double Re_w(const Slice& slice) const {
	const double De = hydraulicDiameter(slice);
	return std::pow(De, 2) * angular_vel() * slice._gas_density / slice._gas_visc;
      }
      
      /*! \brief Heat transfer coefficient between the bed and the gas. 
       */
      double h_gas_bed(const Slice& slice) const {
	const double bed_frac = bedArea(slice) / kilnArea();
	const double Nu = simcem::HT::kiln::Nu_gas_bed_Li_etal_2005(Re_g(slice), Re_w(slice), bed_frac);
	const double De = hydraulicDiameter(slice);
	return Nu * slice._gas_k / De;
      }

      /*! \brief \f$Q_{g\to s}^{cv}\f$. */
      double Q_gas_bed_cv(const Slice& slice) const {
	return h_gas_bed(slice) * chordLength(slice) * (slice._gas->T() - slice._solid->T());
      }
      
      /*! \brief Heat transfer coefficient between the bed and the gas. 
       */
      double h_gas_wall(const Slice& slice) const {
	const double Nu = simcem::HT::kiln::Nu_gas_wall_Li_etal_2005(Re_g(slice), Re_w(slice));
	const double De = hydraulicDiameter(slice);
	return Nu * slice._gas_k / De;
      }

      /*! \brief \f$Q_{g\to s}^{cv}\f$. */
      double Q_gas_wall_cv(const Slice& slice) const {
	return h_gas_wall(slice) * exposedWallPerimeter(slice) * (slice._gas->T() - slice._T_wall);
      }
      
      double h_ext_amb_conv(const Slice& slice) const {
	double amb_visc, amb_k;
	std::tie(amb_visc, amb_k) = simcem::trans::NASA_transport(_ambient);
	const double amb_dens = _ambient->M() / _ambient->V();
	
	const double outerDiam = 2 * outerRadius();
	const double Pr = (_ambient->Cp() / _ambient->M()) * amb_visc / amb_k;
	const double Tf = (_ambient->T() + slice._T_ext_shell) / 2.0;
	const double Gr = simcem::Database::g * std::pow(amb_dens, 2) * (1/Tf) * std::abs(slice._T_ext_shell - _ambient->T()) * std::pow(outerDiam, 3) / std::pow(amb_visc, 2);
	  
	const double Nu = simcem::HT::natConv::Nu_horiz_cyl(Gr, Pr);

	return Nu * amb_k / outerDiam;
      }

      /*! \brief \f$Q_{g\to s}^{cv}\f$. */
      double R_sh_ext_cv(const Slice& slice) const {
	return 1.0 / (h_ext_amb_conv(slice) * outerRadius() * 2 * simcem::Database::pi);
      }

      
      /*! \brief External radiative heat transfer coefficient. */
      double h_ext_amb_rad(const Slice& slice) const {
	return simcem::HT::rad::h_rad_simple(_ambient->T(), slice._T_ext_shell,
					     _shell_emissivity,
					     _ambient->db()->steffanBoltzmannConstant);
      }
      
      double R_sh_ext_rd(const Slice& slice) const {
	return 1.0 / (h_ext_amb_rad(slice) * outerRadius() * 2 * simcem::Database::pi);
      }

      double Q_w_ext(const Slice& slice) const {
	const double Rtotal = R_wall_shell(slice) + 1.0 / (1.0 / R_sh_ext_rd(slice) + 1.0 / R_sh_ext_cv(slice));

	return (slice._T_wall - _ambient->T()) / Rtotal;
      }

      double Q_w_sh(const Slice& slice) const {
	return (slice._T_wall - slice._T_ext_shell) / R_wall_shell(slice);
      }
      
      /*! \brief One-zone radiative network for a kiln slice.

	The wall and bed are considered as two thermally homogeneous
	gray bodies which exchange radiative heat through a gray gas
	(gray implying that we are ignoring the frequency aspects of
	radiation). 

	The same network as Fig.6 in Li_etal_2005, but this is also
	Fig.8-39 in JP Holman's Heat Transfer, 10th Ed.
	
	\image html kiln_radiation_network.svg
	
	The electrical resistance analogy is then used. The junctions
	(circles) are annotated with their "potentials"
	(\f$E_{bed},\,J_b,\,J_w,\,E_{gas}\f$) whereas the other
	annotations are resistances.  As the two nodes \f$J_b\f$ and
	\f$J_w\f$ are merely "virtual" junctions they cannot
	accumulate heat/current, so their heat/current must sum to
	zero:

	\f[
	\frac{E_{bed}-J_b}{\left(1-\varepsilon_b\right)\varepsilon_b^{-1}\,A_b^{-1}} + \frac{J_w-J_b}{\tau_g^{-1}\,F_{bw}^{-1}\,A_b^{-1}} + \frac{E_{gas}-J_b}{\varepsilon_g^{-1}\,F_{bg}^{-1}\,A_b^{-1}}=0
	\f]\f[
	\frac{E_{wall}-J_w}{\left(1-\varepsilon_w\right)\varepsilon_w^{-1}\,A_w^{-1}} + \frac{J_b-J_w}{\tau_g^{-1}\,F_{bw}^{-1}\,A_b^{-1}} + \frac{E_{gas}-J_w}{\varepsilon_g^{-1}\,F_{wg}^{-1}\,A_w^{-1}}=0
	\f]

	The emission fluxes (potentials) are already known:
	\f[
	E_{bed} = \sigma\,T^4_{bed} \qquad E_{gas} = \sigma\,T^4_{gas} \qquad E_{wall} = \sigma\,T^4_{wall}
	\f]

	Thus the task is to determine the radiosities \f$J_{bed}\f$
	and \f$J_{wall}\f$ by solving the two simultaneous equations
	above. Once this is complete the heat fluxes/currents from the wall and
	bed are obtainable from the network's potentials and
	resistances:
	
	\f[
	q_{bed} = \frac{E_{bed}-J_b}{\left(1-\varepsilon_b\right)^{-1}\varepsilon_b\,A_b}
	\qquad
	q_{wall} = \frac{E_{wall}-J_w}{\left(1-\varepsilon_w\right)^{-1}\varepsilon_w\,A_w}
	\f]
	\f[
	q_{gas} = \frac{J_b-E_{gas}}{\varepsilon_g^{-1}\,F_{bg}^{-1}\,A_b^{-1}} + \frac{J_w-E_{gas}}{\varepsilon_g^{-1}\,F_{wg}^{-1}\,A_w^{-1}} = -q_{bed} -q_{wall}
	\f]
	Although the reciprocity relationship
	(\f$F_{12}A_{1}=F_{21}A_{2}\f$) and summation rule (\f$\sum_b
	F_{ab}=1\f$) allow us to change the view factors, the
	selection above has been made as \f$F_{bw}=F_{bg}=F_{wg}=1\f$
	thus these can be eliminated.  For simplicity it is assumed
	that spectral effects on the gas absorptivity (resulting from
	the different temperatures of the wall and gas) can be ignored
	and \f$\tau_g=1-\varepsilon_g\f$. The areas are given in terms
	of per unit length of kiln, thus
	\f$A_w=R\left(2\,\pi-\theta\right)\f$ and \f$A_b\f$ is the
	chord length (Li_etal_2005 incorrectly uses the covered
	perimeter of the inner kiln surface). All the resistance terms
	can then be determined a priori.
	
      */
      std::array<double, 3> rad_fluxes(const Slice& slice) const {
	const double Aw = exposedWallPerimeter(slice);
	const double Abed = chordLength(slice);
	const double h = _bedModel->bedHeight(slice);
	const double pathLength = simcem::HT::kiln::pathLength(_innerRadius, h);
	const Model& gas = *slice._gas;
	const Model& solid = *slice._solid;
	const double e_gas = simcem::HT::rad::Hottel_emissivity(gas, pathLength);
	const double sigma = gas.db()->steffanBoltzmannConstant;

	const std::array<double, 3> V = {{
	    sigma * std::pow(solid.T(), 4),
	    sigma * std::pow(slice._T_wall, 4),
	    sigma * std::pow(gas.T(), 4)
	  }};
	
	const std::array<double, 5> Rinv = {{
	    _bed_emissivity * Abed / (1 - _bed_emissivity),
	    (1 - e_gas) * Abed,
	    e_gas * Abed,
	    _wall_emissivity * Aw / (1 - _wall_emissivity),
	    e_gas * Aw
	  }};
	
	return simcem::HT::rad::twobodygas_network(V, Rinv);
      }	

      void printFluxes(const Slice& slice) const {
	auto rad = rad_fluxes(slice);
	const double Q_s_rad = rad[0], Q_w_rad = rad[1], Q_g_rad = rad[2];
	
	std::cout << "Z=" << slice._Z << ", Tg=" << slice._gas->T() << ", Ts=" << slice._solid->T() << ", Tw=" << slice._T_wall << ", Tshell=" << slice._T_ext_shell << std::endl;
	std::cout << "Q_s_rad=" << Q_s_rad << ", Q_w_rad="<<Q_w_rad<<", Q_g_rad="<<Q_g_rad<< std::endl;
	std::cout << "Q_gs_cv=" << Q_gas_bed_cv(slice) << std::endl;
	std::cout << "Q_ws_cd=" << Q_wall_bed_cd(slice) << std::endl;
	std::cout << "Q_gw_cv=" << Q_gas_wall_cv(slice) << std::endl;
	std::cout << "Q_w_sh=" << Q_w_sh(slice) << std::endl;
	std::cout << "Q_w_ext=" << Q_w_ext(slice) << std::endl;
      }

      
      /*! \brief Solves the system as a steady-state inert system.

	The differential equations to be solved are:
	
	\f[
	\dot{m}_s\,c_{p,s} \frac{\partial T_s}{\partial z} = Q_{g\to s}^{cv}+Q_{w\to s}^{cd}+Q_{g\to s}^{rd}+Q_{w\to s}^{rd}
	\f]
	\f[
	\dot{m}_g\,c_{p,g} \frac{\partial T_g}{\partial z} = -\left(Q_{g\to s}^{cv}+Q_{g\to w}^{cv}+Q_{g\to s}^{rd}+Q_{g\to w}^{rd}\right)
	\f]

	These equations require a solution for the wall \f$T_w\f$ and
	shell \f$T_{shell}\f$ temperatures, which are solved from the
	following algebraic equations:

	\f[
	Q_{sh\to ext} = Q_{w\to sh}
	\f]
	\f[
	Q_{w\to sh} = Q_{g\to w}^{rd} + Q_{g\to w}^{cv} - Q_{w\to s}^{rd} - Q_{w\to s}^{cd}
	\f]
      */
      void solve_SS_inert(Slice slice_in, std::vector<realtype> stop_points, bool store_intermediate) {
	using namespace sundials;
	
	//For debugging, output all kiln settings
	//std::cout << "Kiln setup" << std::endl;
	//std::cout << _RPM << " "
	//	  << _innerRadius << " "
	//	  << _length << " "
	//	  << _particle_diam << " "
	//	  << _shell_emissivity << " "
	//	  << _bed_emissivity << " "
	//	  << _wall_emissivity << " "
	//	  << _solid_density << " "
	//	  << _solid_k << " "
	//	  << _bedModel << " "
	//	  << _ambient << " "
	//	  << std::endl;
	//
	//for (const auto& layer : _layers)
	//  std::cout << "Layer: " << layer._thickness << " " << layer._k << std::endl;
	//
	// TODO: add slice data output
	
	struct SolverData {
	  Kiln& kiln;
	  Slice slice;

	  void save_state(Serial_Vector& T) const {
	    T[0] = slice._gas->T();
	    T[1] = slice._solid->T();
	    T[2] = slice._T_wall;
	    T[3] = slice._T_ext_shell;
	  }

	  void load_state(Serial_Vector& T, const double Z) {
	    slice._gas->set  (simcem::Objective_t::T, T[0], simcem::Objective_t::p);
	    slice._solid->set(simcem::Objective_t::T, T[1], simcem::Objective_t::p);
	    slice._T_wall = T[2];
	    slice._T_ext_shell = T[3];
	    slice.compute_properties(kiln.gasArea(slice), kiln.solid_k(T[1]));
	    slice._Z = Z;
	  }

	  static int resrob(realtype Z, N_Vector T_, N_Vector Tprime_, N_Vector residuals_, void *user_data) {
	    Serial_Vector T(T_), Tprime(Tprime_), residuals(residuals_);
	    SolverData& data = *((SolverData*)user_data);
	    
	    //Load the slice state
	    data.load_state(T, Z);

	    auto rad = data.kiln.rad_fluxes(data.slice);
	    const double Q_s_rad = rad[0], Q_w_rad = rad[1], Q_g_rad = rad[2];

	    //The solid energy balance
	    //-ve as we're starting at the hot end
	    //The gas energy balance
	    residuals[0] = -(-data.kiln.Q_gas_bed_cv(data.slice) - data.kiln.Q_gas_wall_cv(data.slice) + Q_g_rad) - data.slice._gas->Cp() * Tprime[0];

	    residuals[1] = (data.kiln.Q_gas_bed_cv(data.slice) + data.kiln.Q_wall_bed_cd(data.slice) + Q_s_rad) - data.slice._solid->Cp() * Tprime[1];
	    
	    //The shell-wall relationship
	    residuals[2] = Q_w_rad + data.kiln.Q_gas_wall_cv(data.slice) - data.kiln.Q_wall_bed_cd(data.slice) - data.kiln.Q_w_ext(data.slice);
	    residuals[3] = data.kiln.Q_w_sh(data.slice) - data.kiln.Q_w_ext(data.slice);

	    return 0;
	  }
	};
	
	//Set up the system data;
	SolverData data{*this, slice_in};
	
	Serial_Vector T(4);
	data.save_state(T);

	Serial_Vector Tprime(4);
	Tprime[0] = 0;
	Tprime[1] = 0;
	Tprime[2] = Tprime[3] = 0;

	IDA solver(slice_in._Z, SolverData::resrob, T, Tprime);
	
	solver.setTol(1e-4, 1e-2);
	solver.setUserData(data);
	solver.denseSolver(T);

	Serial_Vector Id(4);
	Id[0] = Id[1] = 1.0; //Differential 
	Id[2] = Id[3] = 0.0; //Algebraic
	solver.setId(Id);

	solver.calcIC(IDA_YA_YDP_INIT, _length);

	Serial_Vector T_interp(4);
	Serial_Vector Tprime_interp(4);

	std::sort(stop_points.begin(), stop_points.end());
	auto stop_it = stop_points.begin();
	
	while (stop_it != stop_points.end()) {
	  const realtype Znext = solver.solve(*stop_it, T, Tprime, IDA_ONE_STEP);

	  while (Znext > *stop_it) {
	    solver.interpolate(*stop_it, T_interp, 0);
	    solver.interpolate(*stop_it, Tprime_interp, 1);
	    data.load_state(T_interp, *stop_it);
	    _slices.push_back(data.slice);
	    if (++stop_it == stop_points.end()) break;
	  }

	  if (store_intermediate && (stop_it != stop_points.end())) {
	    data.load_state(T, Znext);
	    _slices.push_back(data.slice);
	  }
	}
      }

      std::vector<Slice>& getSlices() { return _slices; }

    protected:
      
      const double _RPM;
      const double _innerRadius;
      const double _length;
      const double _particle_diam;
      const double _shell_emissivity;
      const double _bed_emissivity;
      const double _wall_emissivity;
      const double _solid_density;
      const double _bed_void_frac;

      sym::Expr _solid_k;
      
      std::vector<Layer> _layers;
      simcem::shared_ptr<BedModel> _bedModel;
      simcem::shared_ptr<simcem::ModelIdealGasTp> _ambient;      
      std::vector<Slice> _slices;
    };
  }
}
