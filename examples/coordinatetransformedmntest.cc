// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <examples/coordinate-transformed-mn.hh>
#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/sourcebeam.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/checkerboard.hh>

int main(int argc, char* argv[])
{
  if (argc > 1)
    DXTC_CONFIG.read_command_line(argc, argv);

  using GridType1d = YaspGrid<1, EquidistantOffsetCoordinates<double, 1>>;
  using GridType3d = YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>;
  using GV1d = typename GridType1d::LeafGridView;
  using GV3d = typename GridType3d::LeafGridView;
  using HFBasis1d = HatFunctionMomentBasis<double, 1, double, 50>;
  using PMBasis1d = PartialMomentBasis<double, 1, double, 50>;
  using MBasis1d = LegendreMomentBasis<double, double, 50>;
  using HFBasis3d = HatFunctionMomentBasis<double, 3, double, /*refinements = */ 2>;
  using PMBasis3d = PartialMomentBasis<double, 3, double, /*refinements = */ 1>;
  using MBasis3d = RealSphericalHarmonicsMomentBasis<double, double, /*order = */ 3, /*fluxdim = */ 3>;
  using SourceBeamSolver1dHFMN = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, HFBasis1d>>;
  using SourceBeamSolver1dPMMN = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, PMBasis1d>>;
  using SourceBeamSolver1dMN = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, MBasis1d>>;
  using CheckerboardSolver3dHFMN = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, HFBasis3d>>;
  using CheckerboardSolver3dPMMN = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, PMBasis3d>>;
  using CheckerboardSolver3dMN = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, MBasis3d>>;
  SourceBeamSolver1dHFMN solver1;
  SourceBeamSolver1dPMMN solver2;
  SourceBeamSolver1dMN solver3;
  CheckerboardSolver3dHFMN solver4;
  CheckerboardSolver3dPMMN solver5;
  CheckerboardSolver3dMN solver6;
  solver1.solve();
  solver2.solve();
  solver3.solve();
  solver4.solve();
  solver5.solve();
  solver6.solve();
} // ... main(...)
