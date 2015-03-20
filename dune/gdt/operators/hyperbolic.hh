// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_HYPERBOLIC_HH
#define DUNE_GDT_OPERATORS_HYPERBOLIC_HH

#include <type_traits>

#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/localevaluation/laxfriedrichs.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forward
template< class AnalyticalFluxImp, class LocalizableFunctionImp,class FVSpaceImp >
class HyperbolicLaxFriedrichs;


namespace internal {

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class FVSpaceImp >
class HyperbolicLaxFriedrichsTraits
{
  static_assert(std::is_base_of< Stuff::GlobalFunctionInterface< typename AnalyticalFluxImp::EntityType,
                                                                 typename FVSpaceImp::RangeFieldType,
                                                                 FVSpaceImp::dimRange,
                                                                 typename FVSpaceImp::DomainFieldType,
                                                                 FVSpaceImp::dimDomain >,
                AnalyticalFluxImp >::value,
                "AnalyticalFluxType has to be derived from Stuff::GlobalFunctionInterface!");
  static_assert(Stuff::is_localizable_function< LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_space< FVSpaceImp >::value,    "FVSpaceImp has to be derived from SpaceInterface!");
public:
  typedef HyperbolicLaxFriedrichs< AnalyticalFluxImp, LocalizableFunctionImp, FVSpaceImp >      derived_type;
  typedef FVSpaceImp                                                                            FVSpaceType;
  typedef typename FVSpaceType::RangeFieldType                                                  RangeFieldType;
  typedef typename FVSpaceType::GridViewType                                                    GridViewType;
  typedef AnalyticalFluxImp                                                                     AnalyticalFluxType;
  typedef LocalizableFunctionImp                                                                LocalizableFunctionType;
  typedef DiscreteFunction< FVSpaceType, Dune::Stuff::LA::CommonDenseVector< RangeFieldType > > FVFunctionType;
  typedef typename FVSpaceType::DomainFieldType                                                 FieldType;
}; // class HyperbolicLaxFriedrichsTraits

} // namespace internal


template< class AnalyticalFluxImp, class LocalizableFunctionImp, class FVSpaceImp >
class HyperbolicLaxFriedrichs
  : public Dune::GDT::OperatorInterface< internal::HyperbolicLaxFriedrichsTraits< AnalyticalFluxImp,
                                                                                  LocalizableFunctionImp,
                                                                                  FVSpaceImp > >
  , public SystemAssembler< FVSpaceImp >
{
  typedef Dune::GDT::OperatorInterface< internal::HyperbolicLaxFriedrichsTraits< AnalyticalFluxImp,
                                                                                 LocalizableFunctionImp,
                                                                                 FVSpaceImp > > OperatorBaseType;
  typedef SystemAssembler< FVSpaceImp >     AssemblerBaseType;

  typedef typename Dune::GDT::LocalEvaluation::LaxFriedrichs::Inner< LocalizableFunctionImp >   NumericalFluxType;
  typedef typename Dune::GDT::LocalOperator::Codim1FV< NumericalFluxType >                      LocalOperatorType;
  typedef typename LocalAssembler::Codim1CouplingFV< LocalOperatorType >                        InnerAssemblerType;

public:
  typedef internal::HyperbolicLaxFriedrichsTraits< AnalyticalFluxImp, LocalizableFunctionImp, FVSpaceImp > Traits;
  typedef typename Traits::FVSpaceType             FVSpaceType;
  typedef typename Traits::GridViewType            GridViewType;
  typedef typename Traits::AnalyticalFluxType      AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::FVFunctionType          FVFunctionType;

  HyperbolicLaxFriedrichs(const std::shared_ptr< const AnalyticalFluxType >& analytical_flux,
                 const LocalizableFunctionType& ratio_dt_dx,
                 const FVSpaceType& fv_space)
    : OperatorBaseType()
    , AssemblerBaseType(fv_space)
    , local_operator_(*analytical_flux, ratio_dt_dx)
    , inner_assembler_(local_operator_)
//    , update_function_(fv_space, "solution")
  {}

using AssemblerBaseType::add;
using AssemblerBaseType::assemble;

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range)
  {
    this->add(inner_assembler_, source, range);
    this->assemble();
  }

private:
  const LocalOperatorType local_operator_;
  const InnerAssemblerType inner_assembler_;
//  FVFunctionType update_function_;
}; // class HyperbolicLaxFriedrichs


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_HYPERBOLIC_HH
