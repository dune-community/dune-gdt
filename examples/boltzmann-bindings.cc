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

#include "boltzmann.hh"

using BoltzmannSolver2d = BoltzmannSolver<2>;
using BoltzmannSolver3d = BoltzmannSolver<3>;

using namespace boost::python;

#include <dune/xt/common/disable_warnings.hh>
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads2d, BoltzmannSolver2d::init, 0, 10)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(next_n_timesteps_overloads2d, BoltzmannSolver2d::next_n_timesteps, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_rhs_overloads2d, BoltzmannSolver2d::apply_rhs_operator, 3, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads3d, BoltzmannSolver3d::init, 0, 10)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(next_n_timesteps_overloads3d, BoltzmannSolver3d::next_n_timesteps, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_rhs_overloads3d, BoltzmannSolver3d::apply_rhs_operator, 3, 6)
#include <dune/xt/common/reenable_warnings.hh>


BOOST_PYTHON_MODULE(boltzmann)
{
  using VectorType = typename BoltzmannSolver3d::VectorType;
  using RangeFieldType = typename BoltzmannSolver3d::RangeFieldType;

  // 2d
  VectorType (BoltzmannSolver2d::*apply_rhs_without_params2d)(VectorType, const double) const =
      &BoltzmannSolver2d::apply_rhs_operator;
  VectorType (BoltzmannSolver2d::*apply_rhs_with_params2d)(VectorType,
                                                           const double,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType) =
      &BoltzmannSolver2d::apply_rhs_operator;

  class_<BoltzmannSolver2d>("BoltzmannSolver2d",
                            init<optional<const std::string,
                                          const size_t,
                                          const size_t,
                                          const bool,
                                          const bool,
                                          const double,
                                          const double,
                                          const double,
                                          const double,
                                          const bool>>())
      .def("init", &BoltzmannSolver2d::init, init_overloads2d())
      .def("linear", &BoltzmannSolver2d::linear)
      .def("solve", &BoltzmannSolver2d::solve)
      .def("next_n_timesteps", &BoltzmannSolver2d::next_n_timesteps, next_n_timesteps_overloads2d())
      .def("reset", &BoltzmannSolver2d::reset)
      .def("finished", &BoltzmannSolver2d::finished)
      .def("apply_kinetic_operator", &BoltzmannSolver2d::apply_kinetic_operator)
      .def("apply_restricted_kinetic_operator", &BoltzmannSolver2d::apply_restricted_kinetic_operator)
      .def("prepare_restricted_operator", &BoltzmannSolver2d::prepare_restricted_operator)
      .def("source_dofs", &BoltzmannSolver2d::restricted_op_input_dofs)
      .def("len_source_dofs", &BoltzmannSolver2d::restricted_op_input_dofs_size)
      .def("apply_rhs_operator", apply_rhs_without_params2d)
      .def("apply_rhs_operator", apply_rhs_with_params2d, apply_rhs_overloads2d())
      .def("set_rhs_operator_parameters", &BoltzmannSolver2d::set_rhs_operator_parameters)
      .def("get_initial_values", &BoltzmannSolver2d::get_initial_values)
      .def("current_time", &BoltzmannSolver2d::current_time)
      .def("set_current_time", &BoltzmannSolver2d::set_current_time)
      .def("set_current_solution", &BoltzmannSolver2d::set_current_solution)
      .def("time_step_length", &BoltzmannSolver2d::time_step_length)
      .def("t_end", &BoltzmannSolver2d::t_end);

  // 3d
  VectorType (BoltzmannSolver3d::*apply_rhs_without_params3d)(VectorType, const double) const =
      &BoltzmannSolver3d::apply_rhs_operator;
  VectorType (BoltzmannSolver3d::*apply_rhs_with_params3d)(VectorType,
                                                           const double,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType) =
      &BoltzmannSolver3d::apply_rhs_operator;

  class_<BoltzmannSolver3d>("BoltzmannSolver3d",
                            init<optional<const std::string,
                                          const size_t,
                                          const size_t,
                                          const bool,
                                          const bool,
                                          const double,
                                          const double,
                                          const double,
                                          const double,
                                          const bool>>())
      .def("init", &BoltzmannSolver3d::init, init_overloads3d())
      .def("linear", &BoltzmannSolver3d::linear)
      .def("solve", &BoltzmannSolver3d::solve)
      .def("next_n_time_steps", &BoltzmannSolver3d::next_n_timesteps, next_n_timesteps_overloads3d())
      .def("reset", &BoltzmannSolver3d::reset)
      .def("finished", &BoltzmannSolver3d::finished)
      .def("apply_kinetic_operator", &BoltzmannSolver3d::apply_kinetic_operator)
      .def("apply_restricted_kinetic_operator", &BoltzmannSolver3d::apply_restricted_kinetic_operator)
      .def("prepare_restricted_operator", &BoltzmannSolver3d::prepare_restricted_operator)
      .def("source_dofs", &BoltzmannSolver3d::restricted_op_input_dofs)
      .def("len_source_dofs", &BoltzmannSolver3d::restricted_op_input_dofs_size)
      .def("apply_rhs_operator", apply_rhs_without_params3d)
      .def("apply_rhs_operator", apply_rhs_with_params3d, apply_rhs_overloads3d())
      .def("set_rhs_operator_parameters", &BoltzmannSolver3d::set_rhs_operator_parameters)
      .def("get_initial_values", &BoltzmannSolver3d::get_initial_values)
      .def("current_time", &BoltzmannSolver3d::current_time)
      .def("set_current_time", &BoltzmannSolver3d::set_current_time)
      .def("set_current_solution", &BoltzmannSolver3d::set_current_solution)
      .def("time_step_length", &BoltzmannSolver3d::time_step_length)
      .def("t_end", &BoltzmannSolver3d::t_end);
}