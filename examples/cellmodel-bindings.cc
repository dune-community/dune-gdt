// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2016)

#include "config.h"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <dune/gdt/test/cellmodel/cellmodel.hh>

using namespace boost::python;
using namespace Dune;

#include <dune/xt/common/disable_warnings.hh>
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(visualize_overloads, CellModelSolver::visualize, 3, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(visualize_pfield_overloads, CellModelSolver::visualize_pfield, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(visualize_ofield_overloads, CellModelSolver::visualize_ofield, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(visualize_stokes_overloads, CellModelSolver::visualize_stokes, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(prepare_pfield_operator_overloads,
                                       CellModelSolver::prepare_pfield_operator,
                                       2,
                                       3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(prepare_ofield_operator_overloads,
                                       CellModelSolver::prepare_ofield_operator,
                                       2,
                                       3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(set_pfield_jacobian_state_overloads,
                                       CellModelSolver::set_pfield_jacobian_state,
                                       2,
                                       3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(set_ofield_jacobian_state_overloads,
                                       CellModelSolver::set_ofield_jacobian_state,
                                       2,
                                       3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_pfield_operator_overloads, CellModelSolver::apply_pfield_operator, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_ofield_operator_overloads, CellModelSolver::apply_ofield_operator, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_stokes_operator_overloads, CellModelSolver::apply_stokes_operator, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_pfield_jacobian_overloads, CellModelSolver::apply_pfield_jacobian, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_ofield_jacobian_overloads, CellModelSolver::apply_ofield_jacobian, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_stokes_jacobian_overloads, CellModelSolver::apply_stokes_jacobian, 1, 2)
#include <dune/xt/common/reenable_warnings.hh>


BOOST_PYTHON_MODULE(cellmodel)
{

  // Boost.Python cannot handle more than about 14 arguments in the constructor by default
  // TODO: Find a way to increase that limit in boost or change constructor.
  class_<CellModelSolver>("CellModelSolver",
                          init<optional<const std::string,
                                        const double,
                                        const unsigned int,
                                        const unsigned int,
                                        const bool,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        // const double,
                                        // const double,
                                        // const double,
                                        const double>>())
      .def("visualize", &CellModelSolver::visualize, visualize_overloads())
      .def("visualize_pfield", &CellModelSolver::visualize_pfield, visualize_pfield_overloads())
      .def("visualize_ofield", &CellModelSolver::visualize_ofield, visualize_ofield_overloads())
      .def("visualize_stokes", &CellModelSolver::visualize_stokes, visualize_stokes_overloads())
      .def("prepare_pfield_operator", &CellModelSolver::prepare_pfield_operator, prepare_pfield_operator_overloads())
      .def("prepare_ofield_operator", &CellModelSolver::prepare_ofield_operator, prepare_ofield_operator_overloads())
      .def("prepare_stokes_operator", &CellModelSolver::prepare_stokes_operator)
      .def("apply_pfield_operator", &CellModelSolver::apply_pfield_operator, apply_pfield_operator_overloads())
      .def("apply_ofield_operator", &CellModelSolver::apply_ofield_operator, apply_ofield_operator_overloads())
      .def("apply_stokes_operator", &CellModelSolver::apply_stokes_operator, apply_stokes_operator_overloads())
      .def("set_pfield_vec", &CellModelSolver::set_pfield_vec)
      .def("set_ofield_vec", &CellModelSolver::set_ofield_vec)
      .def("set_stokes_vec", &CellModelSolver::set_stokes_vec)
      .def("solve", &CellModelSolver::solve)
      .def("next_n_timesteps", &CellModelSolver::next_n_timesteps)
      .def("apply_inverse_pfield_operator", &CellModelSolver::apply_inverse_pfield_operator)
      .def("apply_inverse_ofield_operator", &CellModelSolver::apply_inverse_ofield_operator)
      .def("apply_inverse_stokes_operator", &CellModelSolver::apply_inverse_stokes_operator)
      .def("apply_pfield_product_operator", &CellModelSolver::apply_pfield_product_operator)
      .def("apply_ofield_product_operator", &CellModelSolver::apply_ofield_product_operator)
      .def("apply_stokes_product_operator", &CellModelSolver::apply_stokes_product_operator)
      .def("set_pfield_jacobian_state",
           &CellModelSolver::set_pfield_jacobian_state,
           set_pfield_jacobian_state_overloads())
      .def("set_ofield_jacobian_state",
           &CellModelSolver::set_ofield_jacobian_state,
           set_ofield_jacobian_state_overloads())
      .def("apply_pfield_jacobian", &CellModelSolver::apply_pfield_jacobian, apply_pfield_jacobian_overloads())
      .def("apply_ofield_jacobian", &CellModelSolver::apply_ofield_jacobian, apply_ofield_jacobian_overloads())
      .def("apply_stokes_jacobian", &CellModelSolver::apply_stokes_jacobian, apply_stokes_jacobian_overloads())
      .def("apply_inverse_pfield_jacobian", &CellModelSolver::apply_inverse_pfield_jacobian)
      .def("apply_inverse_ofield_jacobian", &CellModelSolver::apply_inverse_ofield_jacobian)
      .def("apply_inverse_stokes_jacobian", &CellModelSolver::apply_inverse_stokes_jacobian)
      .def("update_pfield_parameters", &CellModelSolver::update_pfield_parameters)
      .def("update_ofield_parameters", &CellModelSolver::update_ofield_parameters)
      .def("num_cells", &CellModelSolver::num_cells)
      .def("finished", &CellModelSolver::finished)
      .def("linear", &CellModelSolver::linear)
      .def("pfield_vec",
           &CellModelSolver::pfield_vec,
           boost::python::return_value_policy<boost::python::return_by_value>())
      .def("ofield_vec",
           &CellModelSolver::ofield_vec,
           boost::python::return_value_policy<boost::python::return_by_value>())
      .def("stokes_vec",
           &CellModelSolver::stokes_vec,
           boost::python::return_value_policy<boost::python::return_by_value>());
}
