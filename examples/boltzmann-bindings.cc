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

#include "boltzmann.hh"

using namespace boost::python;

using VectorType = typename XT::LA::Container<double, XT::LA::Backends::common_dense>::VectorType;

#define BOLTZMANNBINDINGS(SolverName)                                                                                  \
  VectorType (SolverName::*apply_rhs_without_params##SolverName)(VectorType, const double) const =                     \
      &SolverName::apply_rhs_operator;                                                                                 \
  VectorType (SolverName::*apply_rhs_with_params##SolverName)(VectorType, const double, const std::vector<double>&) =  \
      &SolverName::apply_rhs_operator;                                                                                 \
  class_<SolverName>(#SolverName,                                                                                      \
                     init<optional<const std::string,                                                                  \
                                   const size_t,                                                                       \
                                   const size_t,                                                                       \
                                   const bool,                                                                         \
                                   const bool,                                                                         \
                                   const std::vector<double>&,                                                         \
                                   const double>>())                                                                   \
      .def("linear", &SolverName::linear)                                                                              \
      .def("solve", &SolverName::solve)                                                                                \
      .def("next_n_timesteps", &SolverName::next_n_timesteps)                                                          \
      .def("reset", &SolverName::reset)                                                                                \
      .def("finished", &SolverName::finished)                                                                          \
      .def("apply_kinetic_operator", &SolverName::apply_kinetic_operator)                                              \
      .def("apply_restricted_kinetic_operator", &SolverName::apply_restricted_kinetic_operator)                        \
      .def("prepare_restricted_operator", &SolverName::prepare_restricted_operator)                                    \
      .def("source_dofs", &SolverName::restricted_op_input_dofs)                                                       \
      .def("len_source_dofs", &SolverName::restricted_op_input_dofs_size)                                              \
      .def("apply_rhs_operator", apply_rhs_without_params##SolverName)                                                 \
      .def("apply_rhs_operator", apply_rhs_with_params##SolverName)                                                    \
      .def("set_parameters", &SolverName::set_parameters)                                                              \
      .def("get_initial_values", &SolverName::get_initial_values)                                                      \
      .def("current_time", &SolverName::current_time)                                                                  \
      .def("set_current_time", &SolverName::set_current_time)                                                          \
      .def("set_current_solution", &SolverName::set_current_solution)                                                  \
      .def("time_step_length", &SolverName::time_step_length)                                                          \
      .def("t_end", &SolverName::t_end)                                                                                \
      .def("visualize", &SolverName::visualize);


BOOST_PYTHON_MODULE(boltzmann)
{
  using namespace boost::python;
  using namespace Dune;
  using namespace Dune::GDT;
  using GridType1d = YaspGrid<1, EquidistantOffsetCoordinates<double, 1>>;
  using GridType2d = YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>;
  using GridType3d = YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>;
  using GV1d = typename GridType1d::LeafGridView;
  using GV2d = typename GridType2d::LeafGridView;
  using GV3d = typename GridType3d::LeafGridView;
  using HF50Basis1d = HatFunctionMomentBasis<double, 1, double, 50>;
  using PM10Basis1d = PartialMomentBasis<double, 1, double, 10>;
  using PM20Basis1d = PartialMomentBasis<double, 1, double, 20>;
  using PM30Basis1d = PartialMomentBasis<double, 1, double, 30>;
  using PM40Basis1d = PartialMomentBasis<double, 1, double, 40>;
  using PM50Basis1d = PartialMomentBasis<double, 1, double, 50>;
  using M50Basis1d = LegendreMomentBasis<double, double, 50>;
  using HF66Basis3d = HatFunctionMomentBasis<double, 3, double, /*refinements = */ 2>;
  using PM128Basis3d = PartialMomentBasis<double, 3, double, /*refinements = */ 1>;
  using M3Basis3d = RealSphericalHarmonicsMomentBasis<double, double, /*order = */ 3, /*fluxdim = */ 3>;
  using HFM50SourceBeamSolver = BoltzmannSolver<SourceBeamMn<GV1d, HF50Basis1d>>;
  using PMM50SourceBeamSolver = BoltzmannSolver<SourceBeamMn<GV1d, PM50Basis1d>>;
  using M50SourceBeamSolver = BoltzmannSolver<SourceBeamMn<GV1d, M50Basis1d>>;
  using HFM66CheckerboardSolver3d = BoltzmannSolver<CheckerboardMn<GV3d, HF66Basis3d>>;
  using PMM128CheckerboardSolver3d = BoltzmannSolver<CheckerboardMn<GV3d, PM128Basis3d>>;
  using M3CheckerboardSolver3d = BoltzmannSolver<CheckerboardMn<GV3d, M3Basis3d>>;
  using HFM50PlaneSourceSolver = BoltzmannSolver<PlaneSourceMn<GV1d, HF50Basis1d>>;
  using PMM10PlaneSourceSolver = BoltzmannSolver<PlaneSourceMn<GV1d, PM10Basis1d>>;
  using PMM20PlaneSourceSolver = BoltzmannSolver<PlaneSourceMn<GV1d, PM20Basis1d>>;
  using PMM30PlaneSourceSolver = BoltzmannSolver<PlaneSourceMn<GV1d, PM30Basis1d>>;
  using PMM40PlaneSourceSolver = BoltzmannSolver<PlaneSourceMn<GV1d, PM40Basis1d>>;
  using PMM50PlaneSourceSolver = BoltzmannSolver<PlaneSourceMn<GV1d, PM50Basis1d>>;
  using M50PlaneSourceSolver = BoltzmannSolver<PlaneSourceMn<GV1d, M50Basis1d>>;

  BOLTZMANNBINDINGS(HFM50SourceBeamSolver);
  BOLTZMANNBINDINGS(PMM50SourceBeamSolver);
  BOLTZMANNBINDINGS(M50SourceBeamSolver);
  BOLTZMANNBINDINGS(HFM50PlaneSourceSolver);
  BOLTZMANNBINDINGS(PMM10PlaneSourceSolver);
  BOLTZMANNBINDINGS(PMM20PlaneSourceSolver);
  BOLTZMANNBINDINGS(PMM30PlaneSourceSolver);
  BOLTZMANNBINDINGS(PMM40PlaneSourceSolver);
  BOLTZMANNBINDINGS(PMM50PlaneSourceSolver);
  BOLTZMANNBINDINGS(M50PlaneSourceSolver);
  BOLTZMANNBINDINGS(HFM66CheckerboardSolver3d);
  BOLTZMANNBINDINGS(PMM128CheckerboardSolver3d);
  BOLTZMANNBINDINGS(M3CheckerboardSolver3d);
}
