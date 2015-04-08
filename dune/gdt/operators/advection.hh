// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_HH
#define DUNE_GDT_OPERATORS_ADVECTION_HH

#include <type_traits>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/container/interfaces.hh>

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
template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class RangeImp >
class AdvectionLaxFriedrichsLocalizable;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class FVSpaceImp >
class AdvectionLaxFriedrichs;


namespace internal {

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class RangeImp >
class AdvectionLaxFriedrichsLocalizableTraits
{
  static_assert(std::is_base_of< Stuff::GlobalFunctionInterface< typename AnalyticalFluxImp::EntityType,
                                                                 typename SourceImp::RangeFieldType,
                                                                 SourceImp::dimRange,
                                                                 typename SourceImp::DomainFieldType,
                                                                 SourceImp::dimDomain >,
                AnalyticalFluxImp >::value,
                "AnalyticalFluxImp has to be derived from Stuff::GlobalFunctionInterface!");
  static_assert(Stuff::is_localizable_function< LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_discrete_function< SourceImp >::value, "SourceImp has to be derived from DiscreteFunction!");
  static_assert(is_discrete_function< RangeImp >::value, "RangeImp has to be derived from DiscreteFunction!");
public:
  typedef AdvectionLaxFriedrichsLocalizable< AnalyticalFluxImp,
                                             LocalizableFunctionImp,
                                             SourceImp,
                                             RangeImp >                                         derived_type;
  typedef AnalyticalFluxImp                                                                     AnalyticalFluxType;
  typedef LocalizableFunctionImp                                                                LocalizableFunctionType;
  typedef SourceImp                                                                             SourceType;
  typedef RangeImp                                                                              RangeType;
  typedef typename RangeType::SpaceType::GridViewType                                           GridViewType;
  typedef typename GridViewType::ctype                                                          FieldType;
}; // class AdvectionLaxFriedrichsLocalizableTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class FVSpaceImp >
class AdvectionLaxFriedrichsTraits
{
  static_assert(std::is_base_of< Stuff::GlobalFunctionInterface< typename AnalyticalFluxImp::EntityType,
                                                                 typename FVSpaceImp::RangeFieldType,
                                                                 FVSpaceImp::dimRange,
                                                                 typename FVSpaceImp::DomainFieldType,
                                                                 FVSpaceImp::dimDomain >,
                AnalyticalFluxImp >::value,
                "AnalyticalFluxImp has to be derived from Stuff::GlobalFunctionInterface!");
  static_assert(Stuff::is_localizable_function< LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_space< FVSpaceImp >::value,    "FVSpaceImp has to be derived from SpaceInterface!");
public:
  typedef AdvectionLaxFriedrichs< AnalyticalFluxImp, LocalizableFunctionImp, FVSpaceImp >       derived_type;
  typedef AnalyticalFluxImp                                                                     AnalyticalFluxType;
  typedef LocalizableFunctionImp                                                                LocalizableFunctionType;
  typedef FVSpaceImp                                                                            FVSpaceType;
  typedef typename FVSpaceType::GridViewType                                                    GridViewType;
  typedef typename FVSpaceType::DomainFieldType                                                 FieldType;
}; // class LaxFriedrichsTraits

} // namespace internal


template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class RangeImp >
class AdvectionLaxFriedrichsLocalizable
  : public Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp,
                                                                                LocalizableFunctionImp,
                                                                                SourceImp,
                                                                                RangeImp > >
  , public SystemAssembler< typename RangeImp::SpaceType >
{
  typedef Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp,
                                                                                LocalizableFunctionImp,
                                                                                SourceImp,
                                                                                RangeImp > >    OperatorBaseType;
  typedef SystemAssembler< typename RangeImp::SpaceType >                                       AssemblerBaseType;
public:
  typedef internal::AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp,
                                                             LocalizableFunctionImp,
                                                             SourceImp,
                                                             RangeImp >                         Traits;

  typedef typename Traits::GridViewType                                                         GridViewType;
  typedef typename Traits::SourceType                                                           SourceType;
  typedef typename Traits::RangeType                                                            RangeType;
  typedef typename Traits::AnalyticalFluxType                                                   AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType                                              LocalizableFunctionType;

  typedef typename Dune::GDT::LocalEvaluation::LaxFriedrichs::Inner< LocalizableFunctionImp >   NumericalFluxType;
  typedef typename Dune::GDT::LocalEvaluation::LaxFriedrichs::Absorbing< LocalizableFunctionImp >   NumericalBoundaryFluxType;
  typedef typename Dune::GDT::LocalOperator::Codim1FV< NumericalFluxType >                      LocalOperatorType;
  typedef typename Dune::GDT::LocalOperator::Codim1FVBoundary< NumericalBoundaryFluxType >      LocalBoundaryOperatorType;
  typedef typename LocalAssembler::Codim1CouplingFV< LocalOperatorType >                        InnerAssemblerType;
  typedef typename LocalAssembler::Codim1BoundaryFV< LocalBoundaryOperatorType >                BoundaryAssemblerType;

  AdvectionLaxFriedrichsLocalizable(const AnalyticalFluxType& analytical_flux,
                                    const LocalizableFunctionType& ratio_dt_dx,
                                    const SourceType& source,
                                    RangeType& range)
    : OperatorBaseType()
    , AssemblerBaseType(range.space())
    , local_operator_(analytical_flux, ratio_dt_dx)
    , local_boundary_operator_(analytical_flux)
    , inner_assembler_(local_operator_)
    , boundary_assembler_(local_boundary_operator_)
    , source_(source)
    , range_(range)
  {}

  const GridViewType& grid_view() const
  {
    return range_.space().grid_view();
  }

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  RangeType& range()
  {
    return range_;
  }

using AssemblerBaseType::add;
using AssemblerBaseType::assemble;

  void apply()
  {
    this->add(inner_assembler_, source_, range_, new DSG::ApplyOn::InnerIntersections< GridViewType >());
    this->add(inner_assembler_, source_, range_, new DSG::ApplyOn::PeriodicIntersections< GridViewType >());
    this->add(boundary_assembler_, source_, range_, new DSG::ApplyOn::NonPeriodicBoundaryIntersections< GridViewType >());
    this->assemble();
  }

private:
  const LocalOperatorType local_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
  const InnerAssemblerType inner_assembler_;
  const BoundaryAssemblerType boundary_assembler_;
  const SourceType& source_;
  RangeType& range_;
}; // class AdvectionLaxFriedrichsLocalizable



template< class AnalyticalFluxImp, class LocalizableFunctionImp, class FVSpaceImp >
class AdvectionLaxFriedrichs
  : public Dune::GDT::OperatorInterface< internal::AdvectionLaxFriedrichsTraits<  AnalyticalFluxImp,
                                                                                  LocalizableFunctionImp,
                                                                                  FVSpaceImp > >
{
  typedef Dune::GDT::OperatorInterface< internal::AdvectionLaxFriedrichsTraits<  AnalyticalFluxImp,
                                                                                 LocalizableFunctionImp,
                                                                                 FVSpaceImp > > OperatorBaseType;

public:
  typedef internal::AdvectionLaxFriedrichsTraits< AnalyticalFluxImp, LocalizableFunctionImp, FVSpaceImp > Traits;
  typedef typename Traits::GridViewType            GridViewType;
  typedef typename Traits::AnalyticalFluxType      AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::FVSpaceType             FVSpaceType;

  AdvectionLaxFriedrichs(const AnalyticalFluxType& analytical_flux,
                         const LocalizableFunctionType& ratio_dt_dx,
                         const FVSpaceType& fv_space)
    : OperatorBaseType()
    , analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
    , fv_space_(fv_space)
  {}

  const GridViewType& grid_view() const
  {
    return fv_space_.grid_view();
  }

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range) const
  {
    AdvectionLaxFriedrichsLocalizable< AnalyticalFluxType,
                                       LocalizableFunctionType,
                                       SourceType,
                                       RangeType > localizable_operator(analytical_flux_, ratio_dt_dx_, source, range);
    localizable_operator.apply();
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
  const FVSpaceType& fv_space_;
}; // class AdvectionLaxFriedrichs


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_HH
