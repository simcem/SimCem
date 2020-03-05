#include <simcem/simcem.hpp>
#include <simcem/kiln.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
namespace py = pybind11;
using namespace pybind11::literals;

typedef std::vector<simcem::kiln::Slice> SliceList;
PYBIND11_MAKE_OPAQUE(SliceList);
PYBIND11_MAKE_OPAQUE(std::vector<double>);

PYBIND11_MODULE(kiln, m)
{
  py::bind_vector<SliceList>(m, "SliceList");
  py::bind_vector<std::vector<double>>(m, "DoubleList");
  
  py::class_<simcem::kiln::Slice>(m, "Slice")
    .def(py::init<std::shared_ptr<simcem::ModelIdealGasTp>, std::shared_ptr<simcem::ModelIncompressible>, double>(), py::arg("gas"), py::arg("solid"), py::arg("Z0"))
    .def(py::init<const simcem::Components, const double, const double, const double, const double, const double, const double, const double, const std::shared_ptr<simcem::Database> >(), py::arg("solid"), py::arg("volAirFlow"), py::arg("volGasFlow"), py::arg("volO2Flow"), py::arg("volSO2Flow"), py::arg("Tsolid"), py::arg("Tgas"), py::arg("Z0"), py::arg("db"))
    .def_readwrite("gas", &simcem::kiln::Slice::_gas)
    .def_readwrite("solid", &simcem::kiln::Slice::_solid)
    .def_readwrite("T_wall", &simcem::kiln::Slice::_T_wall)
    .def_readwrite("T_ext_shell", &simcem::kiln::Slice::_T_ext_shell)
    .def_readwrite("Z", &simcem::kiln::Slice::_Z)
    ;

  py::class_<simcem::kiln::Kiln>(m, "Kiln")
    .def(py::init<double,double,double,double,double,double,double,double,double, sym::Expr, std::shared_ptr<simcem::Database>>(), py::arg("RPM"), py::arg("innerRadius"), py::arg("length"), py::arg("particleDiam"), py::arg("shell_emissivity"), py::arg("bed_emissivity"), py::arg("wall_emissivity"), py::arg("solid_density"), py::arg("bed_void_frac"), py::arg("solid_k"), py::arg("db"))
    .def("add_layer", &simcem::kiln::Kiln::add_layer, py::arg("material"), py::arg("thickness"), py::arg("k"))
    .def("solve_SS_inert", &simcem::kiln::Kiln::solve_SS_inert, py::arg("init_slice"), py::arg("stop_points"), py::arg("store_intermediate"))
    .def("setFixedHeightBedModel", &simcem::kiln::Kiln::setFixedHeightBedModel, py::arg("solidFraction"))
    .def("getSlices", &simcem::kiln::Kiln::getSlices, py::return_value_policy::reference)
    .def("printFluxes", &simcem::kiln::Kiln::printFluxes)
    .def("R_sh_ext_cv", &simcem::kiln::Kiln::R_sh_ext_cv)
    .def("R_sh_ext_rd", &simcem::kiln::Kiln::R_sh_ext_rd)
    .def("R_wall_shell", &simcem::kiln::Kiln::R_wall_shell)
    .def("Q_w_ext", &simcem::kiln::Kiln::Q_w_ext)
    .def("Q_w_sh", &simcem::kiln::Kiln::Q_w_sh)
    ;
}
