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

#include "coordinate-transformed-mn.hh"

using CoordinateTransformedBoltzmannSolver3d = CoordinateTransformedBoltzmannSolver<3>;

using namespace boost::python;

#include <dune/xt/common/disable_warnings.hh>
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads3d, CoordinateTransformedBoltzmannSolver3d::init, 0, 9)
#include <dune/xt/common/reenable_warnings.hh>

BOOST_PYTHON_MODULE(coordinatetransformedmn)
{
  // 3d
  class_<CoordinateTransformedBoltzmannSolver3d>("CoordinateTransformedBoltzmannSolver3d",
                                                 init<optional<const std::string,
                                                               const size_t,
                                                               const size_t,
                                                               const bool,
                                                               const bool,
                                                               const double,
                                                               const double,
                                                               const double,
                                                               const double>>())
      .def("init", &CoordinateTransformedBoltzmannSolver3d::init, init_overloads3d())
      .def("linear", &CoordinateTransformedBoltzmannSolver3d::linear)
      .def("solve", &CoordinateTransformedBoltzmannSolver3d::solve)
      //   .def("next_n_timesteps", &CoordinateTransformedBoltzmannSolver3d::next_n_timesteps,
      //   next_n_timesteps_overloads3d())
      .def("reset", &CoordinateTransformedBoltzmannSolver3d::reset)
      .def("finished", &CoordinateTransformedBoltzmannSolver3d::finished)
      .def("apply_operator", &CoordinateTransformedBoltzmannSolver3d::apply_operator)
      .def("apply_restricted_operator", &CoordinateTransformedBoltzmannSolver3d::apply_restricted_operator)
      .def("prepare_restricted_operator", &CoordinateTransformedBoltzmannSolver3d::prepare_restricted_operator)
      //   .def("source_dofs", &CoordinateTransformedBoltzmannSolver3d::restricted_op_input_dofs)
      //   .def("len_source_dofs", &CoordinateTransformedBoltzmannSolver3d::restricted_op_input_dofs_size)
      .def("set_parameters", &CoordinateTransformedBoltzmannSolver3d::set_parameters)
      .def("get_initial_values", &CoordinateTransformedBoltzmannSolver3d::get_initial_values)
      .def("current_time", &CoordinateTransformedBoltzmannSolver3d::current_time)
      .def("set_current_time", &CoordinateTransformedBoltzmannSolver3d::set_current_time)
      .def("set_current_solution", &CoordinateTransformedBoltzmannSolver3d::set_current_solution)
      .def("t_end", &CoordinateTransformedBoltzmannSolver3d::t_end);
}
