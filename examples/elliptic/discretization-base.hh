// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_BASE_HH
#define DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_BASE_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/solver/eigen.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localevaluation/sipdg.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/space/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/discretefunction/difference.hh>
#include <dune/gdt/operator/projections.hh>
#include <dune/gdt/operator/products.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>

namespace Discretization {


template< class GridPartType, int polynomialOrder >
class Base
{
public:
  static const unsigned int             dimDomain = GridPartType::dimension;
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int dimRange = 1;
  typedef double            RangeFieldType;
  static const unsigned int polOrder = polynomialOrder;

  typedef Dune::Stuff::FunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::IntersectionType > BoundaryInfoType;

  typedef Dune::Stuff::LA::EigenRowMajorSparseMatrix< RangeFieldType >  MatrixType;
  typedef Dune::Stuff::LA::EigenDenseVector< RangeFieldType >           VectorType;

  Base(const GridPartType& gp,
       const BoundaryInfoType& info,
       const FunctionType& diff,
       const FunctionType& forc,
       const FunctionType& dir,
       const FunctionType& neu)
    : grid_part_(gp)
    , boundary_info_(info)
    , diffusion_(diff)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neu)
  {}

  const GridPartType& grid_part() const
  {
    return grid_part_;
  }

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  const FunctionType& diffusion() const
  {
    return diffusion_;
  }

  const FunctionType& force() const
  {
    return force_;
  }

  const FunctionType& dirichlet() const
  {
    return dirichlet_;
  }

  const FunctionType& neumann() const
  {
    return neumann_;
  }

  template< class ReferenceSolutionType, class DiscreteSolutionType >
  std::vector< RangeFieldType > compute_errors(const ReferenceSolutionType& referencec_solution,
                                               const DiscreteSolutionType& discrete_solution) const
  {
    using namespace Dune;
    using namespace Dune::GDT;
    // checks
    static_assert(std::is_base_of< Stuff::LocalizableFunction, ReferenceSolutionType >::value,
                  "ReferenceSolutionType has to be derived from Stuff::LocalizableFunction!");
    static_assert(std::is_base_of< Stuff::LocalizableFunction, DiscreteSolutionType >::value,
                  "ReferenceSolutionType has to be derived from Stuff::LocalizableFunction!");
    // difference
    typedef DiscreteFunction::Difference< ReferenceSolutionType, DiscreteSolutionType > DifferenceType;
    const DifferenceType difference(referencec_solution, discrete_solution);
    // L2 error
    ProductOperator::L2< GridPartType > l2_product_operator(grid_part_);
    const RangeFieldType l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    // H1 error
    ProductOperator::H1< GridPartType > h1_product_operator(grid_part_);
    const RangeFieldType h1_error = std::sqrt(h1_product_operator.apply2(difference, difference));
    return { l2_error, h1_error };
  }

private:
  const GridPartType& grid_part_;
  const BoundaryInfoType& boundary_info_;
  const FunctionType& diffusion_;
  const FunctionType& force_;
  const FunctionType& dirichlet_;
  const FunctionType& neumann_;
}; // class Base


} // namespace Discretization

#endif // DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_BASE_HH
