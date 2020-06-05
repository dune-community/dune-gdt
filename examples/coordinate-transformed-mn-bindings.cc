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

#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/sourcebeam.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/checkerboard.hh>

#include "coordinate-transformed-mn.hh"

#define COORDINATETRANSFORMEDMNBINDINGS(SolverName)                                                                    \
  class_<SolverName>(#SolverName,                                                                                      \
                     init<optional<const std::string,                                                                  \
                                   const size_t,                                                                       \
                                   const size_t,                                                                       \
                                   const bool,                                                                         \
                                   const bool,                                                                         \
                                   const std::vector<double>>>())                                                      \
      .def("apply_operator", &SolverName::apply_operator)                                                              \
      .def("apply_restricted_operator", &SolverName::apply_restricted_operator)                                        \
      .def("current_time", &SolverName::current_time)                                                                  \
      .def("dx", &SolverName::dx)                                                                                      \
      .def("finished", &SolverName::finished)                                                                          \
      .def("get_initial_values", &SolverName::get_initial_values)                                                      \
      .def("len_source_dofs", &SolverName::restricted_op_input_dofs_size)                                              \
      .def("linear", &SolverName::linear)                                                                              \
      .def("prepare_restricted_operator", &SolverName::prepare_restricted_operator)                                    \
      .def("set_current_solution", &SolverName::set_current_solution)                                                  \
      .def("set_current_time", &SolverName::set_current_time)                                                          \
      .def("set_parameters", &SolverName::set_parameters)                                                              \
      .def("solve", &SolverName::solve)                                                                                \
      .def("source_dofs", &SolverName::restricted_op_input_dofs)                                                       \
      .def("t_end", &SolverName::t_end)                                                                                \
      .def("u_from_alpha", &SolverName::u_from_alpha)                                                                  \
      .def("visualize", &SolverName::visualize)

BOOST_PYTHON_MODULE(coordinatetransformedmn)
{
  using namespace boost::python;
  using namespace Dune;
  using namespace Dune::GDT;
  using GridType1d = YaspGrid<1, EquidistantOffsetCoordinates<double, 1>>;
  using GridType3d = YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>;
  using GV1d = typename GridType1d::LeafGridView;
  using GV3d = typename GridType3d::LeafGridView;
  using HF50Basis1d = HatFunctionMomentBasis<double, 1, double, 50>;
  using PM50Basis1d = PartialMomentBasis<double, 1, double, 50>;
  using M50Basis1d = LegendreMomentBasis<double, double, 50>;
  using HF66Basis3d = HatFunctionMomentBasis<double, 3, double, /*refinements = */ 2>;
  using PM128Basis3d = PartialMomentBasis<double, 3, double, /*refinements = */ 1>;
  using M3Basis3d = RealSphericalHarmonicsMomentBasis<double, double, /*order = */ 3, /*fluxdim = */ 3>;
  using HFM50SourceBeamSolver = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, HF50Basis1d>>;
  using PMM50SourceBeamSolver = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, PM50Basis1d>>;
  using M50SourceBeamSolver = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, M50Basis1d>>;
  using HFM66CheckerboardSolver3d = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, HF66Basis3d>>;
  using PMM128CheckerboardSolver3d = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, PM128Basis3d>>;
  using M3CheckerboardSolver3d = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, M3Basis3d>>;

  COORDINATETRANSFORMEDMNBINDINGS(HFM50SourceBeamSolver);
  COORDINATETRANSFORMEDMNBINDINGS(PMM50SourceBeamSolver);
  COORDINATETRANSFORMEDMNBINDINGS(M50SourceBeamSolver);
  COORDINATETRANSFORMEDMNBINDINGS(HFM66CheckerboardSolver3d);
  COORDINATETRANSFORMEDMNBINDINGS(PMM128CheckerboardSolver3d);
  COORDINATETRANSFORMEDMNBINDINGS(M3CheckerboardSolver3d);
}
