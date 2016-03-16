// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_BASE_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_BASE_HH


namespace Dune {
namespace GDT {
namespace Hyperbolic {


enum class ChooseDiscretizer
{
  fv
}; // enum class ChooseDiscretizer

enum class FluxTimeStepperKombinations
{
  godunov_euler,
  godunov_adaptiveRK,
  laxfriedrichs_euler,
  godunovwithreconstruction_euler
};


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_DISCRETIZERS_BASE_HH
