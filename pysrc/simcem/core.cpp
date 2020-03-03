#include "simcem/simcem.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace simcem;

// Get all of the opaque containers defined
typedef std::unordered_map<std::string, simcem::Element> ElementMap;
PYBIND11_MAKE_OPAQUE(ElementMap)

typedef std::unordered_map<std::string, simcem::Component> ComponentMap;
PYBIND11_MAKE_OPAQUE(ComponentMap)

typedef std::vector<simcem::Isotope> IsotopeList;
PYBIND11_MAKE_OPAQUE(IsotopeList)

typedef std::vector<simcem::Component::Property> PropertyList;
PYBIND11_MAKE_OPAQUE(PropertyList)

PYBIND11_MAKE_OPAQUE(simcem::Component::PropertyMap)

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


PYBIND11_MODULE(core, m)
{
  m.def("solveCubic", solveCubic);

  py::bind_map<ElementMap>(m, "ElementMap");
  py::bind_map<ComponentMap>(m, "ComponentMap");
  py::bind_vector<IsotopeList>(m, "IsotopeList");
  py::bind_vector<PropertyList>(m, "PropertyList");
  
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
    .def_readwrite("g", &Database::g)
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

  py::bind_map<Components>(m, "Components")
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

    py::class_<Data>(m, "Data")
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

  py::class_<Component::Property, Data>(m, "Property")
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

}
