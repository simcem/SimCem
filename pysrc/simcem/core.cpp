#include "simcem/simcem.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace simcem;

PYBIND11_MODULE(core, m)
{
  py::class_<simcem::Database, std::shared_ptr<simcem::Database> >(m, "Database")
    .def(py::init(&Database::create))
    .def(py::init(&Database::create_from_file))
    //.def("load", &Database::load)
    //.def("xml", &Database::xml)
    //.def("getComponent", &Database::getComponent, py::return_value_policy::reference)
    //.def("getComponents", &Database::getComponents, py::return_value_policy::reference)
    //.def("getElement", &Database::getElement, py::return_value_policy::reference)
    //.def("getElements", &Database::getElements, py::return_value_policy::reference)
    //.def("registerComponent", &Database::registerComponent, py::return_value_policy::reference)
    //.def_readwrite("avogadro", &Database::avogadro)
    //.def_readwrite("kB", &Database::kB)
    //.def_readwrite("steffanBoltzmannConstant", &Database::steffanBoltzmannConstant)
    //.def_readwrite("g", &Database::g)
    //.def_property_readonly("R", &Database::R)
   ;

}
