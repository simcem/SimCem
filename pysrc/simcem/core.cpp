#include <simcem/simcem.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace simcem;

// Get all of the opaque containers defined
typedef std::unordered_map<std::string, simcem::Element> ElementMap;
PYBIND11_MAKE_OPAQUE(ElementMap);
typedef std::unordered_map<std::string, simcem::Component> ComponentMap;
PYBIND11_MAKE_OPAQUE(ComponentMap);
typedef std::vector<simcem::Isotope> IsotopeList;
PYBIND11_MAKE_OPAQUE(IsotopeList);
typedef std::vector<simcem::Component::Property> PropertyList;
PYBIND11_MAKE_OPAQUE(PropertyList);
PYBIND11_MAKE_OPAQUE(std::vector<FunctionCurve::ShomateTerm>);
PYBIND11_MAKE_OPAQUE(std::vector<Bond>);
PYBIND11_MAKE_OPAQUE(std::vector<Atom>);
PYBIND11_MAKE_OPAQUE(std::vector<TabulatedCurve::Datum>);
PYBIND11_MAKE_OPAQUE(simcem::Component::PropertyMap);

Components init_from_dict(py::dict d) {
  Components retval;
  for (auto item : d)
    if (!py::isinstance<py::str>(item.first))
      throw pybind11::key_error();
    else if (py::isinstance<py::float_>(item.second))
      retval[item.first.cast<std::string>()] = item.second.cast<double>();
    else if (py::isinstance<py::int_>(item.second))
      retval[item.first.cast<std::string>()] = item.second.cast<double>();
    else
      throw pybind11::value_error();
  return retval;
}

py::list solveCubic(double A, double B, double C, double D) {
  using namespace sym;
  Polynomial<1> x{0, 1};
  auto roots = solve_real_roots(expand(A * x*x*x + B * x*x + C * x + D));
  py::list retval;
  for (const auto& root : roots)
    retval.append(root);
  return retval;
}

template<class T, class M>
py::class_<T> registerEnum(M& m, std::string name) {
  auto e = py::class_<T>(m, name.c_str())
    .def("__str__", &T::operator std::string)
    .def("__repr__", &T::operator std::string)
    .def("__int__", &T::operator int)
    .def(py::self == py::self)
    ;
  
  for (size_t i(0); i < T::size(); ++i)
    e.attr(T::strings[i]) = T(i);

  return e;
}


template<class T>
std::string xmlstring(const T& t) {
  std::ostringstream os;
  stator::xml::Document doc;
  stator::xml::Node root = doc.add_node("Root");
  t.xml(root);
  root.firstNode().print(os);
  return os.str();
}

//This provides the overloads for concrete classes that have
//overloaded Model's pure virtual functions
template<class ModelBase>
class PyModel : public ModelBase {
public:
  using ModelBase::ModelBase;

  Objective_t Y(size_t idx) const override { PYBIND11_OVERLOAD(Objective_t, ModelBase, Y, idx); }
  double p() const override { PYBIND11_OVERLOAD(double, ModelBase, p, ); }
  double T() const override { PYBIND11_OVERLOAD(double, ModelBase, T, ); }
  double v(std::string molID) const override { PYBIND11_OVERLOAD(double, ModelBase, v, molID); }
  double V() const override { PYBIND11_OVERLOAD_PURE(double, ModelBase, V, ); }
  double s(std::string molID) const override { PYBIND11_OVERLOAD(double, ModelBase, s, molID); }
  double chemPot(std::string molID) const override { PYBIND11_OVERLOAD(double, ModelBase, chemPot, molID); }
  double h(std::string molID) const override { PYBIND11_OVERLOAD(double, ModelBase, h, molID); }
  double u(std::string molID) const override { PYBIND11_OVERLOAD(double, ModelBase, u, molID); }
  double a(std::string molID) const override { PYBIND11_OVERLOAD(double, ModelBase, a, molID); }
  double Cp() const override { PYBIND11_OVERLOAD(double, ModelBase, Cp, ); }
  double Alpha() const override { PYBIND11_OVERLOAD(double, ModelBase, Alpha, ); }
  double Beta() const override { PYBIND11_OVERLOAD(double, ModelBase, Beta, ); }
  std::string str() const override { PYBIND11_OVERLOAD(std::string, ModelBase, str, ); }
  double dfdT(Objective_t obj) const override { PYBIND11_OVERLOAD(double, ModelBase, dfdT, obj); }
  double dfdY2(Objective_t obj) const override { PYBIND11_OVERLOAD(double, ModelBase, dfdY2, obj); }
  double dfdNi(Objective_t obj, std::string molID) const override { PYBIND11_OVERLOAD(double, ModelBase, dfdNi, obj, molID); }
  void set(Objective_t val, double target, Objective_t constant) override { PYBIND11_OVERLOAD(void, ModelBase, set, val, target, constant); }
};

//This specialisation is for the pure virtual methods of the Model base class
template<>
class PyModel<Model> : public Model {
public:
  using Model::Model;

  Objective_t Y(size_t idx) const override { PYBIND11_OVERLOAD_PURE(Objective_t, Model, Y, idx); }
  double p() const override { PYBIND11_OVERLOAD_PURE(double, Model, p, ); }
  double T() const override { PYBIND11_OVERLOAD_PURE(double, Model, T, ); }
  double v(std::string molID) const override { PYBIND11_OVERLOAD_PURE(double, Model, v, molID); }
  double V() const override { PYBIND11_OVERLOAD_PURE(double, Model, V, ); }
  double s(std::string molID) const override { PYBIND11_OVERLOAD_PURE(double, Model, s, molID); }
  double chemPot(std::string molID) const override { PYBIND11_OVERLOAD_PURE(double, Model, chemPot, molID); }
  double h(std::string molID) const override { PYBIND11_OVERLOAD_PURE(double, Model, h, molID); }
  double u(std::string molID) const override { PYBIND11_OVERLOAD_PURE(double, Model, u, molID); }
  double a(std::string molID) const override { PYBIND11_OVERLOAD_PURE(double, Model, a, molID); }
  double Cp() const override { PYBIND11_OVERLOAD_PURE(double, Model, Cp, ); }
  double Alpha() const override { PYBIND11_OVERLOAD_PURE(double, Model, Alpha, ); }
  double Beta() const override { PYBIND11_OVERLOAD_PURE(double, Model, Beta, ); }
  std::string str() const override { PYBIND11_OVERLOAD_PURE(std::string, Model, str, ); }
  double dfdT(Objective_t obj) const override { PYBIND11_OVERLOAD_PURE(double, Model, dfdT, obj); }
  double dfdY2(Objective_t obj) const override { PYBIND11_OVERLOAD_PURE(double, Model, dfdY2, obj); }
  double dfdNi(Objective_t obj, std::string molID) const override { PYBIND11_OVERLOAD_PURE(double, Model, dfdNi, obj, molID); }
  void set(Objective_t val, double target, Objective_t constant) override { PYBIND11_OVERLOAD_PURE(void, Model, set, val, target, constant); }
};


PYBIND11_MODULE(core, m)
{
  m.def("solveCubic", solveCubic);

  py::bind_map<ElementMap>(m, "ElementMap");
  py::bind_map<ComponentMap>(m, "ComponentMap");
  py::bind_vector<IsotopeList>(m, "IsotopeList");
  py::bind_vector<PropertyList>(m, "PropertyList");
  py::bind_vector<std::vector<FunctionCurve::ShomateTerm>>(m, "ShomateTerms");
  py::bind_vector<std::vector<Bond>>(m, "BondList");
  py::bind_vector<std::vector<Atom>>(m, "AtomList");
  py::bind_vector<std::vector<TabulatedCurve::Datum>>(m, "TabulatedCurveList");
  py::bind_map<simcem::Component::PropertyMap>(m, "PropertyMap");
  
  py::class_<simcem::Database, std::shared_ptr<simcem::Database> >(m, "Database")
    .def(py::init(&Database::create))
    .def(py::init(&Database::create_from_file))
    .def("load", &Database::load)
    .def("xml", &Database::xml)
    .def("getComponent", &Database::getComponent, py::return_value_policy::reference)
    .def("getComponents", &Database::getComponents, py::return_value_policy::reference)
    .def("getElement", &Database::getElement, py::return_value_policy::reference)
    .def("getElements", &Database::getElements, py::return_value_policy::reference)
    .def("registerComponent", &Database::registerComponent, py::return_value_policy::reference)
    .def_readwrite("avogadro", &Database::avogadro)
    .def_readwrite("kB", &Database::kB)
    .def_readwrite("steffanBoltzmannConstant", &Database::steffanBoltzmannConstant)
    .def_property_readonly_static("g", +[](){ return Database::g; })
    .def_property_readonly("R", &Database::R)
    ;

  py::class_<simcem::Isotope>(m, "Isotope")
    .def(py::init<size_t, size_t, double, double, double, std::string, std::string, std::string>(), py::arg("Z"), py::arg("N"), py::arg("mass"), py::arg("mass_uncertainty"), py::arg("abundance"), py::arg("symbol"), py::arg("name"), py::arg("category"))
    .def_readonly("Z", &Isotope::_Z)
    .def_readonly("N", &Isotope::_N)
    .def_readonly("mass", &Isotope::_mass)
    .def_readonly("mass_uncertainty", &Isotope::_mass_uncertainty)
    .def_readonly("abundance", &Isotope::_abundance)
    .def_readonly("symbol", &Isotope::_symbol)
    .def_readonly("name", &Isotope::_name)
    .def_readonly("category", &Isotope::_category)
    ;
  
  py::class_<simcem::Element, simcem::Isotope, std::vector<Isotope> >(m, "Element")
    .def(py::init<size_t, size_t, double, double, double, std::string, std::string, std::string, IsotopeList, size_t, size_t, std::string, std::string>(), py::arg("Z"), py::arg("N"), py::arg("mass"), py::arg("mass_uncertainty"), py::arg("abundance"), py::arg("symbol"), py::arg("name"), py::arg("category"), py::arg("isotopes"), py::arg("group"), py::arg("period"), py::arg("block")="", py::arg("referenceComponentID")="")
    .def_readwrite("referenceComponentID", &Element::_referenceComponentID)
    .def_readwrite("group", &Element::_group)
    .def_readwrite("period", &Element::_period)
    .def_readwrite("block", &Element::_block)
    ;

  py::bind_map<Components, shared_ptr<Components>>(m, "Components")
    .def(py::init(&init_from_dict))
    .def("N", &Components::N)
    .def("M", &Components::M)
    .def("m", &Components::m)
    .def("removeSmallComponents", &Components::removeSmallComponents)
    .def("elements", &Components::elements)
    .def("ref_elements", &Components::ref_elements)
    .def("normalised", &Components::normalised)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def(py::self * double())
    .def(double() * py::self)
    .def(py::self / double())
    .def(py::self == py::self)
    .def(py::self != py::self)
    ;

  m.def("MassToMoles", &Components::fromMasses);

  py::class_<sym::Expr>(m, "Expr")
    .def(py::init<std::string>())
    .def(py::init<int>())
    .def(py::init<double>())
    .def("__repr__", +[](const sym::Expr& self) { return stator::repr(self); })
    .def("__str__", +[](const sym::Expr& self) { return stator::repr(self); })
    .def("latex", +[](const sym::Expr& self) { return stator::repr<stator::ReprConfig<stator::Latex_output> >(self); })
    .def("simplify", +[](const sym::Expr& self) { return sym::simplify(self); })
    .def("__add__", +[](const sym::Expr& l, const sym::Expr& r) { return sym::Expr(l+r); })
    .def("__radd__", +[](const sym::Expr& l, const sym::Expr& r) { return sym::Expr(l+r); })
    .def("__sub__", +[](const sym::Expr& l, const sym::Expr& r) { return sym::Expr(l-r); })
    .def("__rsub__", +[](const sym::Expr& l, const sym::Expr& r) { return sym::Expr(l-r); })
    .def("__mul__", +[](const sym::Expr& l, const sym::Expr& r) { return sym::Expr(l*r); })
    .def("__rmul__", +[](const sym::Expr& l, const sym::Expr& r) { return sym::Expr(l*r); })
    .def("__div__", +[](const sym::Expr& l, const sym::Expr& r) { return sym::Expr(l/r); })
    .def("__rdiv__", +[](const sym::Expr& l, const sym::Expr& r) { return sym::Expr(l/r); })
    .def(py::self / py::self)
    ;

  py::implicitly_convertible<int, sym::Expr>();
  py::implicitly_convertible<double, sym::Expr>();

  {
    double (Data::*normalise_units_overload1)(Objective_t variable, double value) const = &Data::normalise_units;
    double (Data::*restore_units_overload1)(Objective_t variable, double value) const = &Data::restore_units;
    std::string (Data::*LaTeX_units_overload1)(Objective_t variable) const = &Data::LaTeX_units;
    double (Data::*normalise_units_overload2)(Property_t variable, double value) const = &Data::normalise_units;
    double (Data::*restore_units_overload2)(Property_t variable, double value) const = &Data::restore_units;
    std::string (Data::*LaTeX_units_overload2)(Property_t variable) const = &Data::LaTeX_units;

    py::class_<Data, shared_ptr<Data>>(m, "Data")
      .def(py::init<>())
      .def("normalise_units", normalise_units_overload1)
      .def("restore_units", restore_units_overload1)
      .def("LaTeX_units", LaTeX_units_overload1)
      .def("normalise_units", normalise_units_overload2)
      .def("restore_units", restore_units_overload2)
      .def("LaTeX_units", LaTeX_units_overload2)
      ;

    py::class_<Data::Comment>(m, "Comment")
      .def_readonly("src", &Data::Comment::_src)
      .def_readonly("text", &Data::Comment::_text)
      ;
  }

  py::class_<Component::Property, Data, std::shared_ptr<Component::Property>>(m, "Property")
    .def(py::init<Property_t, double, std::string>())
    .def_property_readonly("value", &Component::Property::getValue)
    .def_property_readonly("orig_value", &Component::Property::getOrigValue)
    .def_property_readonly("type", &Component::Property::getType)
    .def_property_readonly("source", &Component::Property::getSource)
    ;

  registerEnum<simcem::Objective_t>(m, "Objective_t");
  registerEnum<simcem::Reference_t>(m, "Reference_t");
  registerEnum<simcem::Property_t>(m, "Property_t");
  registerEnum<simcem::T_unit_t>(m, "T_unit_t")
    .def("scale", &simcem::T_unit_t::scale)
    .def("origin", &simcem::T_unit_t::origin)
    ;
  registerEnum<simcem::Q_unit_t>(m, "Q_unit_t")
    .def("scale", &simcem::Q_unit_t::scale)
    ;

  registerEnum<simcem::E_unit_t>(m, "E_unit_t")
    .def("scale", &simcem::E_unit_t::scale)
    ;
  
  registerEnum<simcem::P_unit_t>(m, "P_unit_t")
    .def("scale", &simcem::P_unit_t::scale)
    ;

  registerEnum<simcem::L_unit_t>(m, "L_unit_t")
    .def("scale", &simcem::L_unit_t::scale)
    ;

  py::class_<simcem::Curve, simcem::Data, shared_ptr<simcem::Curve>>(m, "Curve")
    .def("inRange", &Curve::inRange)
    .def("xmin", &Curve::xmin)
    .def("xmax", &Curve::xmax)
    .def_readwrite("comments", &Curve::_comments)
    .def_readwrite("source", &Curve::_source)
    .def_readwrite("reference", &Curve::_reference)
    .def_readwrite("variable", &Curve::_variable)
    .def_readwrite("xvar", &Curve::_xvar)
    .def_readwrite("error", &Curve::_error)
    ;

  py::class_<Isobar, shared_ptr<Isobar>>(m, "Isobar")
    .def(py::init<shared_ptr<Curve>, double>())
    .def("__str__", &xmlstring<Isobar>)
    .def("__repr__", &xmlstring<Isobar>)
    .def_readwrite("security", &Isobar::security)
    .def_property_readonly("p", &Isobar::p)
    .def_property("curve", &Isobar::getCurve, &Isobar::setCurve)
    .def("eval", +[](const shared_ptr<Isobar>& iso, double T){
	auto data = sym::ad<2>(iso->getCurve()->getFunction(), sym::Var<sym::vidx<'T'>>() = T);
	data[2] *= 2;
	return std::vector<double>(&data[0], &data[0]+3);
      })
    ;

  py::class_<FunctionCurve::ShomateTerm>(m, "ShomateTerm")
    .def(py::init<double, double>())
    ;
  
  py::class_<FunctionCurve, Curve, shared_ptr<FunctionCurve> >(m, "FunctionCurve")
    .def(py::init<double, double, std::string, std::string, Objective_t, Reference_t, Objective_t, T_unit_t, Q_unit_t, E_unit_t, P_unit_t, L_unit_t, sym::Expr>())
    .def("LaTeX", &FunctionCurve::LaTeX)
    .def_static("Shomate", &FunctionCurve::Shomate)
    ;

  py::class_<TabulatedCurve::Datum>(m, "TabulatedCurveDatum")
    .def_readonly("x", &TabulatedCurve::Datum::x)
    .def_readonly("val", &TabulatedCurve::Datum::val)
    ;
  
  py::class_<Atom, shared_ptr<Atom> >(m, "Atom")
    .def(py::init<size_t, size_t, int, double, double, double>())
    .def_readwrite("ID", &Atom::_ID)
    .def_readwrite("Z", &Atom::_Z)
    .def_readwrite("N", &Atom::_N)
    .def_readwrite("charge", &Atom::_charge)
    .def_readwrite("bonds", &Atom::_bonds)
    .def_readwrite("x", &Atom::_x)
    .def_readwrite("y", &Atom::_y)
    ;
  
  py::class_<Bond>(m, "Bond")
    .def(py::init<shared_ptr<Atom>, shared_ptr<Atom>, double>())
    .def_readwrite("atom1", &Bond::_atom1)
    .def_readwrite("atom2", &Bond::_atom2)
    .def_readwrite("order", &Bond::_order)
    ;

  const Component::PropertyMap& (Component::*getprop_overload)() const = &Component::getProperties;
  std::vector<shared_ptr<Atom> >&  (Component::*getstruct_overload)() = &Component::getStructure;
    
  py::class_<Component>(m, "Component")
    .def("__str__", &xmlstring<Component>)
    .def("__repr__", &xmlstring<Component>)
    .def("__len__", &Component::size)
    .def("__getitem__", [](const Component &s, size_t i) {
	if (i >= s.size()) throw py::index_error();
	return s[i];
      })
    .def("__iter__", [](const Component &s) { return py::make_iterator(s.begin(), s.end()); },
	 py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    .def("getElements", &Component::getElements, py::return_value_policy::reference_internal)
    .def_property("formula", &Component::getFormula, &Component::setFormula, py::return_value_policy::reference_internal)
    .def("getPhase", &Component::getPhase, py::return_value_policy::reference_internal)
    .def("registerPhase", &Component::registerPhase, py::return_value_policy::reference_internal)
    .def("add_alias", &Component::add_alias)
    .def("remove_alias", &Component::remove_alias)
    .def("getProperties", getprop_overload, py::return_value_policy::reference_internal)
    .def("getStructure", getstruct_overload, py::return_value_policy::reference_internal)
    .def("registerProperty", &Component::registerProperty, py::arg("type"), py::arg("source"), py::arg("val"), py::arg("TUnit")=T_unit_t(T_unit_t::K), py::arg("QUnit")=Q_unit_t(Q_unit_t::mol), py::arg("EUnit")=E_unit_t(E_unit_t::J), py::arg("PUnit")=P_unit_t(P_unit_t::Pa), py::arg("LUnit")=L_unit_t(L_unit_t::m))
    .def("mass", &Component::mass)
    .def("getAliases", &Component::getAliases, py::return_value_policy::reference_internal)
    ;
  
  py::class_<Phase>(m, "Phase")
    .def("__len__", &Phase::size)
    .def("__getitem__", [](const Phase &s, size_t i) {
	if (i >= s.size()) throw py::index_error();
	return s[i];
      })
    .def("__iter__", [](const Phase &s) { return py::make_iterator(s.begin(), s.end()); },
	 py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    .def("registerIsobar", &Phase::registerIsobar)
    .def("add_alias", &Phase::add_alias)
    .def("remove_alias", &Phase::remove_alias)
    .def_readwrite("name", &Phase::_name)
    .def_readonly("type", &Phase::_type)
    .def("getAliases", &Phase::getAliases, py::return_value_policy::reference)
    .def_property("comments", &Phase::getComments, &Phase::setComments)
    ;

  py::class_<Model, Components, PyModel<Model>, shared_ptr<Model>>(m, "Model")
    //Virtual functions
    .def("Y", &Model::Y)
    .def("p", &Model::p)
    .def("T", &Model::T)
    .def("v", &Model::v)
    .def("V", &Model::V)
    .def("s", &Model::s)
    .def("chemPot", &Model::chemPot)
    .def("h", &Model::h)
    .def("u", &Model::u)
    .def("a", &Model::a)
    .def("Cp", &Model::Cp)
    .def("Alpha", &Model::Alpha)
    .def("Beta", &Model::Beta)
    .def("str", &Model::str)
    .def("__repr__", &Model::str)
    .def("__str__", &Model::str)
    .def("dfdT", &Model::dfdT)
    .def("dfdY2", &Model::dfdY2)
    .def("dfdNi", &Model::dfdNi)
    .def("set", &Model::set)
    //Normal functions
    .def("S", &Model::S)
    .def("negs", &Model::s)
    .def("G", &Model::G)
    .def("H", &Model::H)
    .def("U", &Model::U)
    .def("A", &Model::A)
    .def("db", &Model::db)
    .def("f", &Model::f)
    .def("elements", &Model::elements)
    .def("M", &Model::M)
    .def("m", &Model::m)
    .def_property_readonly("components", +[](const Model& m) { return Components(m); })
    .def("hasData", &Model::hasData)
    ;
  
  py::class_<ModelExcess, PyModel<ModelExcess>, Model, shared_ptr<ModelExcess>>(m, "ModelExcess")
    .def(py::init<shared_ptr<Model>>())
    ;

  py::class_<ModelIdealGasTp, PyModel<ModelIdealGasTp>, Model, shared_ptr<ModelIdealGasTp>>(m, "ModelIdealGasTp")
    .def(py::init<shared_ptr<Database>, Components, double, double>(), py::arg("db"), py::arg("components"), py::arg("T")=298.15, py::arg("p")=1e5)
    ;

  py::class_<ModelIdealGasTV, PyModel<ModelIdealGasTV>, Model, shared_ptr<ModelIdealGasTV>>(m, "ModelIdealGasTV")
    .def(py::init<shared_ptr<Database>, Components, double, double>(),py::arg("db"), py::arg("components"), py::arg("T")=298.15, py::arg("V"))
    ;

  py::class_<ModelIncompressible, PyModel<ModelIncompressible>, Model, shared_ptr<ModelIncompressible>>(m, "ModelIncompressible")
    .def(py::init<shared_ptr<Database>, Components, double, double, std::string, bool>(), py::arg("db"), py::arg("components"), py::arg("T")=298.15, py::arg("p")=1e5, py::arg("type")="solid", py::arg("mixing")=false)
    ;

  py::class_<ModelIncompressibleSolid, PyModel<ModelIncompressibleSolid>, Model, shared_ptr<ModelIncompressibleSolid>>(m, "ModelIncompressibleSolid")
    .def(py::init<shared_ptr<Database>, Components, double, double>(), py::arg("db"), py::arg("components"), py::arg("T")=298.15, py::arg("p")=1e5)
    ;

  py::class_<ModelIncompressibleLiquid, PyModel<ModelIncompressibleLiquid>, Model, shared_ptr<ModelIncompressibleLiquid>>(m, "ModelIncompressibleLiquid")
    .def(py::init<shared_ptr<Database>, Components, double, double>(), py::arg("db"), py::arg("components"), py::arg("T")=298.15, py::arg("p")=1e5)
    ;

  
  py::bind_vector<MultiPhase>(m, "MultiPhase")
    .def("T", &MultiPhase::T)
    .def("N", &MultiPhase::T)
    .def("p", &MultiPhase::p)
    .def("M", &MultiPhase::M)
    .def("G", &MultiPhase::G)
    .def("H", &MultiPhase::H)
    .def("U", &MultiPhase::U)
    .def("V", &MultiPhase::V)
    .def("S", &MultiPhase::S)
    .def("A", &MultiPhase::A)
    .def("Cp", &MultiPhase::Cp)
    .def("Alpha", &MultiPhase::Alpha)
    .def("Beta", &MultiPhase::Beta)
    ;


  py::enum_<System::Optimiser_t>(m, "Optimiser_t")
    .value("SLSQP", System::Optimiser_t::SLSQP)
    .value("NelderMead", System::Optimiser_t::NelderMead)
    .value("IPopt", System::Optimiser_t::IPopt)
    ;

  py::class_<System, MultiPhase>(m, "System")
    .def(py::init<Objective_t, Objective_t, bool, System::Optimiser_t, double, size_t, bool, bool, bool>(), py::arg("Y1"), py::arg("Y2"), py::arg("reactive")=false, py::arg("alg")=System::Optimiser_t::IPopt, py::arg("molarConstraintTol")=1e-8, py::arg("max_eval")=100, py::arg("debug")=false, py::arg("scale_problem")=true, py::arg("rank_reduce_constraints")=true)
    .def("__repr__", &System::str)
    .def("__str__", &System::str)
    .def("equilibrate", &System::equilibrate, py::arg("Y1init") =HUGE_VAL, py::arg("Y2init") =HUGE_VAL)
    .def("setEnsemble", &System::setEnsemble)
    ;
  
  m.def("NASA_transport", simcem::trans::NASA_transport);
}
