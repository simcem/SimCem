#pragma once
#include <simcem/enum.hpp>
#include <simcem/config.hpp>
#include <boost/algorithm/string.hpp>
#include <stator/symbolic/symbolic.hpp>
#include <algorithm>

namespace simcem {
  struct Reference_t : Enum<Reference_t> {
    using Enum<Reference_t>::Enum;
    enum { Elements298, Elements273, Elements, Self };
    static const constexpr char* const strings[] = {"Elements298", "Elements273", "Elements", "Self"};
  };

  struct Property_t : Enum<Property_t> {
    using Enum<Property_t>::Enum;
    enum { vc, pc, AcentricFactor, Tc};
    static const constexpr char* const strings[] = {"vc", "pc", "AcentricFactor", "Tc"};
  };
  
  struct Objective_t : Enum<Objective_t> {
    using Enum<Objective_t>::Enum;
    enum { G, H, negS, S, U, A, p, V, T, Cp, Cv, dynvisc, thermcond, density};
    static const constexpr char* const strings[] = {"G", "H", "negS", "S", "U", "A", "p", "V", "T", "Cp", "Cv", "dynvisc", "thermcond", "density"};
  };

  struct T_unit_t : Enum<T_unit_t> {
    using Enum<T_unit_t>::Enum;
    enum { K, C, F, R };
    static const constexpr char* const strings[] = {"K", "C", "F", "R"};
    inline double scale() const {
      static constexpr double data[] = {1.0, 1.0, 5.0/9.0, 5.0/9.0};
      return data[_value];
    }
    inline double origin() const {
      static constexpr double data[] = {0, 273.15, 273.15 - 32.0 * 5.0 / 9.0, 0};
      return data[_value];
    }
  };

  struct Q_unit_t : Enum<Q_unit_t> {
    using Enum<Q_unit_t>::Enum;
    enum { mol, kmol, g, kg, lb};
    static const constexpr char* const strings[] = {"mol", "kmol", "g", "kg", "lb"};
    inline double scale() const {
      static constexpr double data[] = {1.0, 1000.0, 1e-3, 1, 0.45359237};
      return data[_value];
    }
  };

  struct E_unit_t : Enum<E_unit_t> {
    using Enum<E_unit_t>::Enum;
    enum { J, kJ, cal, kcal };
    static const constexpr char* const strings[] = {"J", "kJ", "cal", "kcal"};
    inline double scale() const {
      static constexpr double data[] = {1.0, 1000.0, 4.184, 4184.0};
      return data[_value];
    }
  };

  struct P_unit_t : Enum<P_unit_t> {
    using Enum<P_unit_t>::Enum;
    enum { Pa, bar, at, atm, torr, psi };
    static const constexpr char* const strings[] = {"Pa", "bar", "at", "atm", "torr", "psi"};
    inline double scale() const {
      static constexpr double data[] = {1.0, 1e5, 9.80665e4, 1.01325e5, 101325.0/760, 6.8948e3 };
      return data[_value];
    }
  };

  struct L_unit_t : Enum<L_unit_t> {
    using Enum<L_unit_t>::Enum;
    enum { m, cm, mm, ft, inches};
    static const constexpr char* const strings[] = {"m", "cm", "mm", "ft", "inches"};
    inline double scale() const {
      static constexpr double data[] = {1.0, 0.01, 0.001, 0.3048, 0.3048/12};
      return data[_value];
    }
  };
  
  struct Data {
    Data():
      _T("K"),
      _Q("mol"),
      _E("J"),
      _P("Pa"),
      _L("m"),
      _mass(0)
    {}
    
    void xml(stator::xml::Node node) const {      
      if (!_source.empty())
	node.add_attribute("Source", _source);
      if (_T != T_unit_t::K)
	node.add_attribute("TUnit", std::string(_T));
      if (_Q != Q_unit_t::mol)
	node.add_attribute("QUnit", std::string(_Q));
      if (_E != E_unit_t::J)
	node.add_attribute("EUnit", std::string(_E));
      if (_P != P_unit_t::Pa)
	node.add_attribute("PUnit", std::string(_P));
      if (_L != L_unit_t::m)
	node.add_attribute("LUnit", std::string(_L));

      for (const auto& comment : _comments) {
	//Don't output empty comments
	if (stator::strip(comment._text).empty() &&
	    (stator::strip(comment._src).empty() || (stator::strip(comment._src) == stator::strip(_source))))
	    continue;
	auto comnode = node.add_node("Comment");
	if (!comment._src.empty())
	  comnode.add_attribute("Source", comment._src);
	comnode = comment._text;
      }
    }

    Data(stator::xml::Node node):
      _T("K"),
      _Q("mol"),
      _E("J"),
      _P("Pa"),
      _L("m"),
      _mass(0)
    {
      //Only load values if present
      for (stator::xml::Node comnode = node.findNode("Comment"); comnode.valid(); ++comnode) {
	std::string src;
	if (comnode.hasAttribute("Source"))
	  src = comnode.getAttribute("Source").as<std::string>();
	_comments.emplace_back(src, comnode);
      }
      
      if (node.hasAttribute("Source"))
	_source = node.getAttribute("Source").as<std::string>();

      if (node.hasAttribute("TUnit"))
	_T = T_unit_t(node.getAttribute("TUnit"));
      if (node.hasAttribute("QUnit"))
	_Q = Q_unit_t(node.getAttribute("QUnit"));
      if (node.hasAttribute("EUnit"))
	_E = E_unit_t(node.getAttribute("EUnit"));
      if (node.hasAttribute("PUnit"))
	_P = P_unit_t(node.getAttribute("PUnit"));
      if (node.hasAttribute("LUnit"))
	_L = L_unit_t(node.getAttribute("LUnit"));
    }

    Data(std::string comment, std::string source,
	 T_unit_t Tunit, Q_unit_t Qunit, E_unit_t Eunit, P_unit_t Punit, L_unit_t Lunit):
      _comments({{std::string(""),comment}}),
      _source(source),
      _T(Tunit), _Q(Qunit), _E(Eunit), _P(Punit), _L(Lunit),
      _mass(0)
    {}

    double normalise_units(Objective_t variable, double value) const {
      const double origin = (variable == Objective_t::T) ? _T.origin() : 0.0;
      return value * scale(variable) + origin;
    }

    double restore_units(Objective_t variable, double value) const {
      const double origin = (variable == Objective_t::T) ? _T.origin() : 0.0;
      return (value - origin) / scale(variable);
    }
    
    double scale(Objective_t variable) const {
      double scale = 1;

      double inv_qscale = 1.0 / _Q.scale();
      if ((_Q == Q_unit_t::g) || (_Q == Q_unit_t::kg) || (_Q == Q_unit_t::lb)) {
	if (_mass == 0)
	  stator_throw() << "Cannot use Isobars with mass units without first registering them to a phase.";
	inv_qscale *= _mass;
      }
      
      switch (variable) {
      case Objective_t::G:
      case Objective_t::H:
      case Objective_t::U:
      case Objective_t::A:
	return _E.scale() * inv_qscale;
      case Objective_t::negS:
      case Objective_t::S:
      case Objective_t::Cp:
      case Objective_t::Cv:
	return _E.scale() * inv_qscale / _T.scale();
      case Objective_t::T:
	return _T.scale();
      case Objective_t::p:
	return _P.scale(); 
      case Objective_t::V:
	return std::pow(_L.scale(), 3) * inv_qscale;
      case Objective_t::dynvisc:
	return _P.scale(); //i.e. Pa s, but the time scale is always seconds!
      case Objective_t::thermcond:
	return _E.scale() / _L.scale() / _T.scale();
      case Objective_t::density:
	return 1.0 / (inv_qscale * std::pow(_L.scale(), 3));
      default:
	stator_throw() << "Can't handle units for this yet " + std::string(variable);
      };
      
      return scale;
    }

    std::string LaTeX_units(Objective_t variable) const {
      auto wrapunit = [](std::string u) { return "\\mathrm{"+u+"}"; };
      switch (variable) {
      case Objective_t::G:
      case Objective_t::H:
      case Objective_t::U:
      case Objective_t::A:
	return wrapunit(_E) + "\\," + wrapunit(_Q)+"^{-1}";
      case Objective_t::negS:
      case Objective_t::S:
      case Objective_t::Cp:
      case Objective_t::Cv:
	return wrapunit(_E) + "\\," + wrapunit(_Q)+"^{-1}" + "\\," + wrapunit(_T)+"^{-1}";
      case Objective_t::T:
	return wrapunit(_T);
      case Objective_t::p:
	return wrapunit(_P); 
      case Objective_t::V:
	return wrapunit(_L)+"^3" + "\\," + wrapunit(_Q)+"^{-1}";
      case Objective_t::dynvisc:
	return wrapunit(_P)+"\\,"+wrapunit("s");
      case Objective_t::thermcond:
	return wrapunit(_E)+"\\,"+wrapunit(_L)+"^{-1}\\,"+wrapunit(_T)+"^{-1}";
      case Objective_t::density:
	return wrapunit(_Q)+"\\,"+wrapunit(_L)+"^{-3}";
      default:
	stator_throw() << "Can't handle units for this yet " + std::string(variable);
      }
    }
    
    double normalise_units(Property_t variable, double value) const {
      const double origin = (variable == Property_t::Tc) ? _T.origin() : 0.0;
      return value * scale(variable) + origin;
    }

    double restore_units(Property_t variable, double value) const {
      const double origin = (variable == Property_t::Tc) ? _T.origin() : 0.0;
      return (value - origin) / scale(variable);
    }
    
    double scale(Property_t variable) const {
      switch (variable) {
      case Property_t::vc:
	return std::pow(_L.scale(), 3) / _Q.scale();	
      case Property_t::Tc:
	return _T.scale();
      case Property_t::AcentricFactor:
	return 1.0;
      case Property_t::pc:
	return _P.scale();
      default:
	stator_throw() << "Can't handle units for this yet " + std::string(variable);
      };
    }

    std::string LaTeX_units(Property_t variable) const {
      auto wrapunit = [](std::string u) { return "\\mathrm{"+u+"}"; };
      switch (variable) {
      case Property_t::vc:
	return wrapunit(_L)+"^3\\,"+wrapunit(_Q)+"^{-1}";
      case Property_t::Tc:
	return wrapunit(_T);
      case Property_t::AcentricFactor:
	return "";
      case Property_t::pc:
	return wrapunit(_P);
      default:
	stator_throw() << "Can't handle units for this yet " + std::string(variable);
      }
    }

    void setMass(double mass) { _mass = mass; }

    std::string getSource() { return _source; }

    bool operator==(const Data& d) const {
      return (_comments == d._comments)
	&& (_source == d._source)
	&& (_T == d._T)
	&& (_Q == d._Q)
	&& (_E == d._E)
	&& (_P == d._P)
	&& (_L == d._L)
	&& (_mass == d._mass)
	;
    }

    struct Comment {
      Comment() {}
      
      Comment(std::string src, std::string text):
	_src(src), _text(text)
      {}
      
      std::string _src;
      std::string _text;
      bool operator==(const Comment& b) const {
	return (_src == b._src) && (_text == b._text);
      }
    };
    
    std::vector<Comment> _comments;
    std::string _source;
    T_unit_t _T;
    Q_unit_t _Q;
    E_unit_t _E;
    P_unit_t _P;
    L_unit_t _L;
    double _mass;
  };

  /*! \brief Data structure for a phase of a particular molecule. */
  struct Curve : Data {
    Curve(std::string comment, std::string source, Objective_t variable, Reference_t reference, Objective_t xvar,
	  T_unit_t Tunit, Q_unit_t Qunit, E_unit_t Eunit, P_unit_t Punit, L_unit_t Lunit):
      Data(comment, source, Tunit, Qunit, Eunit, Punit, Lunit),
      _xvar(xvar),
      _variable(variable),
      _reference(reference)
    {}
      
    Curve(Node xml, Objective_t xvar):
      Data(xml),
      _xvar(xvar),
      _variable(xml.getAttribute("Variable")),
      _reference(xml.getAttribute("Reference"))
    {
      if (xml.hasAttribute("Error"))
	_error = xml.getAttribute("Error");
    }

    virtual ~Curve() {}
	  
    void xml(stator::xml::Node node) const {
      node.add_attribute("Variable", std::string(_variable));

      if (!_error.empty())
	node.add_attribute("Error", _error);
      
      Data::xml(node);

      node.add_attribute("Reference", std::string(_reference));
      xml_extended(node);
    }
    
    bool inRange(double x) const {
      return (x >= xmin()) && (x <= xmax());
    }

    virtual double xmin() const = 0;
    virtual double xmax() const = 0;
    
    virtual sym::Expr getFunction() const = 0;
    
    virtual void xml_extended(Node) const = 0;

    Objective_t _xvar;
    Objective_t _variable;
    Reference_t _reference;
    double _mass;
    std::string _error;

    Objective_t getVariable() const { return _variable; }
    
    static shared_ptr<Curve> load(Node, Objective_t);
  };

      
  struct FunctionCurve: public Curve {
    struct ShomateTerm {
      ShomateTerm(double c, double p): C(c), power(p) {}
      bool operator==(const ShomateTerm&c) const { return (c.C == C) && (power == c.power); }
      double C; double power; };
    
    static std::string Shomate(const std::vector<ShomateTerm>& terms, double HConst = 0, double SConst = 0) {
      std::string input;
      
      static auto signwrap = [](std::string a) { return (((a[0] != '-') && (a[0] != '+')) ? "+" : "") + a;  };
      
      for (const auto& t : terms) {
	std::string coeff = signwrap(stator::repr(t.C));
	if (t.power == -1)
	  input = input + coeff + "*(ln(T)+1)";
	else if (t.power == 0)
	  input = input + coeff + "*T*(1-ln(T))";
	else
	  input = input + signwrap(stator::repr(-t.C / (t.power * (t.power + 1))))+"*T^"+stator::repr(t.power+1);
      }

      if (HConst)
	input = input + "+" + stator::repr(HConst);

      if (SConst)
	input = input + "-T*" + stator::repr(SConst);

      return input;
    }
    
    
    static std::string parseOldXML(Node xml) {
      if (xml.getAttribute("Type").as<std::string>() == "Shomate") {
	std::vector<ShomateTerm> terms;
	
	for (Node tnode = xml.findNode("Term"); tnode.valid(); ++tnode) {
	  const std::string type = tnode.getAttribute("Type");
	  if (type.substr(1,1) != "^") stator_throw() << "Shomate only supports polynomial terms";
	  double C = tnode.getAttribute("C").as<double>();
	  double power = boost::lexical_cast<double>(tnode.getAttribute("Type").as<std::string>().substr(2));
	  terms.push_back(ShomateTerm{C, power});
	}

	double HConst = 0, SConst = 0;
	if (xml.hasNode("HConst"))
	  HConst = xml.getNode("HConst").getAttribute("Value").as<double>();
	if (xml.hasNode("SConst"))
	  SConst = xml.getNode("SConst").getAttribute("Value").as<double>();

	return Shomate(terms, HConst, SConst);
      }
      
      static auto signwrap = [](std::string a) { return (((a[0] != '-') && (a[0] != '+')) ? "+" : "") + a;  };
      
      std::string input = "";
      if (xml.hasNode("Function"))
	input = stator::strip(xml.getNode("Function").getValue());
      
      for (Node tnode = xml.findNode("Term"); tnode.valid(); ++tnode) {
	const std::string type = tnode.getAttribute("Type");
	std::string s_C = tnode.getAttribute("C").as<std::string>();
	s_C = signwrap(s_C);
	
	if (type.substr(1,1) == "^") {
	  std::string s_power = tnode.getAttribute("Type").as<std::string>().substr(2);
	  input = input + s_C + "*T^" + s_power;
	} else if (type == "lnx") {
	  input = input + s_C +"*ln(T)";
	} else if (type == "xlnx") {
	  input = input + s_C + "*T*ln(T)";
	} else if (type == "exp") {
	  input = input + s_C + "*exp("+parseOldXML(tnode)+")";
	} else
	  stator_throw() << "Unknown term type";
      }
      
      return input;
    }
    
    FunctionCurve(Node xml, Objective_t xvar):
      Curve(xml, xvar),
      _x_min(xml.getAttribute("xMin").as<double>()),
      _x_max(xml.getAttribute("xMax").as<double>())
    {
      sym::Var<sym::vidx<'T'>> T;
      
      _f = sym::Expr(parseOldXML(xml));
      auto df = sym::simplify(sym::derivative(_f, T));

      if (xml.getAttribute("Type").as<std::string>() == "Shomate") {
	double Tref = 298.15;
	if (xml.hasNode("TRef"))
	  Tref = xml.getNode("TRef").getAttribute("Value").as<double>();
	
	if (xml.hasNode("HRef") && xml.hasNode("SRef")) {
	  auto H_Ref = xml.getNode("HRef").getAttribute("Value").as<double>();
	  auto S_Ref = xml.getNode("SRef").getAttribute("Value").as<double>();

	  auto G_Shift = sym::sub(_f, T = Tref);
	  auto S_Shift = sym::sub(df, T = Tref);

	  _f = _f + H_Ref - G_Shift + Tref * S_Shift - T * (S_Shift + S_Ref);
	}
      }
    }

    virtual double xmin() const {
      return normalise_units(_xvar, _x_min);
    }
    
    virtual double xmax() const {
      return normalise_units(_xvar, _x_max);
    }
	
    FunctionCurve(double xmin, double xmax, std::string comment, std::string source, Objective_t variable, Reference_t reference,
		  Objective_t xvar,
		  T_unit_t Tunit, Q_unit_t Qunit, E_unit_t Eunit, P_unit_t Punit, L_unit_t Lunit, sym::Expr f
		  ):
      Curve(comment, source, variable, reference, xvar, Tunit, Qunit, Eunit, Punit, Lunit),
      _x_min(xmin), _x_max(xmax), _f(f)
    {}

    virtual void xml_extended(Node xml) const {
      xml.add_attribute("Type", "Function");
      xml.add_attribute("xMin", _x_min);
      xml.add_attribute("xMax", _x_max);
      xml.add_node("Function") = stator::repr(_f);
    }

    virtual sym::Expr getFunction() const {
      sym::Var<sym::vidx<'T'>> T;
      //Standardise the units of the function
      return scale(_variable) * sym::sub(_f , T = (T - _T.origin()) * (1 / scale(Objective_t(Objective_t::T))));
    }
    
    virtual std::array<std::string,3> LaTeX() const {
      sym::Expr df = sym::simplify(sym::derivative(_f, sym::VarRT('T')));
      sym::Expr ddf = sym::simplify(sym::derivative(df, sym::VarRT('T')));
    
      return {{stator::repr<stator::ReprConfig<stator::Latex_output> >(_f),
	    stator::repr<stator::ReprConfig<stator::Latex_output> >(df),
	    stator::repr<stator::ReprConfig<stator::Latex_output> >(ddf)
	    }};
    }
    
    double _x_min;
    double _x_max;
    sym::Expr _f;
  };


  struct TabulatedCurve: public Curve {
    TabulatedCurve(std::string comment, std::string source, Objective_t variable, Reference_t reference, Objective_t xvar,
		   T_unit_t Tunit, Q_unit_t Qunit, E_unit_t Eunit, P_unit_t Punit, L_unit_t Lunit):
      Curve(comment, source, variable, reference, xvar, Tunit, Qunit, Eunit, Punit, Lunit)
    {}

    TabulatedCurve(Node xml, Objective_t xvar):
      Curve(xml, xvar)
    {
      std::string input = xml.getNode("Data");
      std::vector<std::string> lines;
      boost::split(lines, input, boost::is_any_of("\n"), boost::token_compress_on);

      for (size_t line(0); line < lines.size(); ++line) {
	std::vector<std::string> values;
	boost::trim_right(lines[line]);
	boost::trim_left(lines[line]);
	boost::split(values, lines[line], boost::is_any_of("\t "), boost::token_compress_on);
	if (lines[line] == "")
	  continue;
	if (values.size() == 0)
	  continue;
	if (values[0][0] == '!')
	  continue;
	if ((values.size() != 2) && (values.size() != 3)) {
	  std::ostringstream os;
	  for (const auto & split: values)
	    os << "\"" << split << "\",";
	  
	  stator_throw() << "Tabulated Curve has bad formatting on line " << line+1 << ", xml path " << xml.getPath() << "\nLine content is :\"" << lines[line] << "\"" << "\n" << "string was split as follows " << os.str()
			 << "\nXML:\n" << xml.print();
	  
	}
	if (values.size() == 2)
	  add_point(boost::lexical_cast<double>(values[0]), boost::lexical_cast<double>(values[1]));
	if (values.size() == 3)
	  add_point(boost::lexical_cast<double>(values[0]), boost::lexical_cast<double>(values[1]), boost::lexical_cast<double>(values[2]));
      }
      if  (_values.empty())
	stator_throw() << "Loaded an empty Tabulated Curve? xml path " << xml.getPath() 
		       << "\nXML:\n" << xml.print();
    }

    virtual double xmin() const {
      if (_values.empty())
	return +HUGE_VAL;
      return normalise_units(_xvar, _values.front().x);
    }
    
    virtual double xmax() const {
      if (_values.empty())
	return -HUGE_VAL;
      return normalise_units(_xvar, _values.back().x);
    }

    virtual void xml_extended(Node xml) const {
      xml.add_attribute("Type", "Tabulated");

      std::string data;
      for (const auto& datum : _values) {      
	data += "\n" + boost::lexical_cast<std::string>(datum.x) + " " + boost::lexical_cast<std::string>(datum.val);
	if (datum.error != 0)
	  data += " " + boost::lexical_cast<std::string>(datum.error);
      }

      xml.add_node("Data") = data;
    }

    virtual sym::Expr getFunction() const {
      stator_throw() << "Not implemented yet";
      //See commented code below for ideas
    }
    
    //virtual DataPoint eval_worker(double x) const {
    //  for (size_t i(0); i < _values.size() - 1; ++i)
    //	if ((x >= _values[i].x) && (x <= _values[i+1].x)) {
    //	  const double dfdx  = (_values[i+1].val - _values[i].val) / (_values[i+1].x - _values[i].x);
    //	  const double val = _values[i].val + dfdx * (x - _values[i].x);
    //	  return DataPoint{x, val, dfdx, 0.0};
    //	}
    //  
    //  if (x == _values.back().x) {
    //	double dfdx = 0;
    //	if (_values.size() - 1)
    //	  dfdx = (_values[_values.size() - 1].val - _values[_values.size() - 2].val)
    //	    / (_values[_values.size() - 1].x - _values[_values.size() - 2].x);
    //	return DataPoint{x, _values.back().val, dfdx, 0.0};
    //  }
    //  
    //  stator_throw() << "Out of range access for Tabulated Curve";
    //}

    void add_point(double x, double val, double error = 0.0) {
      Datum v{x, val, error};
      _values.insert(std::upper_bound(_values.begin(), _values.end(), v), v);
    }
    
    struct Datum {
      bool operator<(const Datum& d) const { return x < d.x; }
      double x;
      double val;
      double error;
    };
    std::vector<Datum> _values;
  };

  inline
  shared_ptr<Curve> Curve::load(Node xml, Objective_t xvar) {
    try {
      std::string type = xml.getAttribute("Type").as<std::string>();
      if ((type == "Function") || (type == "Shomate"))
	return make_shared<FunctionCurve>(xml, xvar);
      else if (type == "Tabulated")
	return make_shared<TabulatedCurve>(xml, xvar);
      else
	stator_throw() << "Unknown Curve type " << type;
    } catch (const stator::Exception& e) {
      stator_throw() << "Failed parsing Curve\n" << xml.getPath() << "\n" << e.what();
    }
  }
}
