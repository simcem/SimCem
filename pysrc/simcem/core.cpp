#include "simcem/simcem.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace simcem;


typedef std::unordered_map<std::string, simcem::Element> ElementMap;
PYBIND11_MAKE_OPAQUE(ElementMap)

typedef std::unordered_map<std::string, simcem::Component> ComponentMap;
PYBIND11_MAKE_OPAQUE(ComponentMap)

typedef std::vector<simcem::Isotope> IsotopeList;
PYBIND11_MAKE_OPAQUE(IsotopeList)

PYBIND11_MODULE(core, m)
{
  py::bind_map<ElementMap>(m, "ElementMap");
  py::bind_map<ComponentMap>(m, "ComponentMap");
  py::bind_vector<IsotopeList>(m, "IsotopeList");
  
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
  
}
