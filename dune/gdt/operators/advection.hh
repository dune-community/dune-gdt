// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_HH
#define DUNE_GDT_OPERATORS_ADVECTION_HH

#include <type_traits>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/localevaluation/godunov.hh>
#include <dune/gdt/localevaluation/laxfriedrichs.hh>
#include <dune/gdt/localevaluation/laxwendroff.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>

#include <dune/gdt/playground/spaces/dg/pdelabproduct.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forwards
template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionLaxFriedrichsLocalizable;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionLaxFriedrichs;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionGodunovLocalizable;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionGodunov;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionLaxWendroffLocalizable;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionLaxWendroff;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionGodunovWithReconstructionLocalizable;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionGodunovWithReconstruction;


namespace internal {

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionLaxFriedrichsLocalizableTraits
{
  static_assert(std::is_base_of< Stuff::GlobalFunctionInterface< typename AnalyticalFluxImp::EntityType,
                                                                 typename SourceImp::RangeFieldType,
                                                                 SourceImp::dimRange,
                                                                 typename SourceImp::DomainFieldType,
                                                                 SourceImp::dimRange,
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
                                             BoundaryValueImp,
                                             RangeImp >                derived_type;
  typedef AnalyticalFluxImp                                            AnalyticalFluxType;
  typedef LocalizableFunctionImp                                       LocalizableFunctionType;
  typedef SourceImp                                                    SourceType;
  typedef BoundaryValueImp                                             BoundaryValueType;
  typedef RangeImp                                                     RangeType;
  typedef typename SourceType::RangeFieldType                          RangeFieldType;
  typedef typename RangeType::SpaceType::GridViewType                  GridViewType;
  typedef typename GridViewType::ctype                                 FieldType;
}; // class AdvectionLaxFriedrichsLocalizableTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionLaxFriedrichsTraits
{
  static_assert(std::is_base_of< Stuff::GlobalFunctionInterface< typename AnalyticalFluxImp::EntityType,
                                                                 typename FVSpaceImp::RangeFieldType,
                                                                 FVSpaceImp::dimRange,
                                                                 typename FVSpaceImp::DomainFieldType,
                                                                 FVSpaceImp::dimRange,
                                                                 FVSpaceImp::dimDomain >,
                                 AnalyticalFluxImp >::value,
                "AnalyticalFluxImp has to be derived from Stuff::GlobalFunctionInterface!");
  static_assert(Stuff::is_localizable_function< LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_space< FVSpaceImp >::value, "FVSpaceImp has to be derived from SpaceInterface!");
public:
  typedef AdvectionLaxFriedrichs< AnalyticalFluxImp, LocalizableFunctionImp,BoundaryValueImp, FVSpaceImp > derived_type;
  typedef AnalyticalFluxImp                                                                     AnalyticalFluxType;
  typedef LocalizableFunctionImp                                                                LocalizableFunctionType;
  typedef BoundaryValueImp                                                                      BoundaryValueType;
  typedef FVSpaceImp                                                                            FVSpaceType;
  typedef typename FVSpaceType::GridViewType                                                    GridViewType;
  typedef typename FVSpaceType::DomainFieldType                                                 FieldType;
}; // class LaxFriedrichsTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionGodunovLocalizableTraits
    : public AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp, LocalizableFunctionImp, SourceImp, BoundaryValueImp, RangeImp >
{
public:
  typedef AdvectionGodunovLocalizable< AnalyticalFluxImp,
                                             LocalizableFunctionImp,
                                             SourceImp,
                                             BoundaryValueImp,
                                             RangeImp >           derived_type;
}; // class AdvectionGodunovLocalizableTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionGodunovTraits
    : public AdvectionLaxFriedrichsTraits< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueImp, FVSpaceImp >
{
public:
  typedef AdvectionGodunov< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueImp, FVSpaceImp > derived_type;
}; // class AdvectionGodunovTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionLaxWendroffLocalizableTraits
    : public AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp, LocalizableFunctionImp, SourceImp, BoundaryValueImp, RangeImp >
{
public:
  typedef AdvectionLaxWendroffLocalizable< AnalyticalFluxImp,
                                             LocalizableFunctionImp,
                                             SourceImp,
                                             BoundaryValueImp,
                                             RangeImp >           derived_type;
}; // class AdvectionLaxWendroffLocalizableTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionLaxWendroffTraits
    : public AdvectionLaxFriedrichsTraits< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueImp, FVSpaceImp >
{
public:
  typedef AdvectionLaxWendroff< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueImp, FVSpaceImp > derived_type;
}; // class AdvectionLaxWendroffTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionGodunovWithReconstructionLocalizableTraits
    : public AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp, LocalizableFunctionImp, SourceImp, BoundaryValueImp, RangeImp >
{
public:
  typedef AdvectionGodunovWithReconstructionLocalizable< AnalyticalFluxImp,
                                             LocalizableFunctionImp,
                                             SourceImp,
                                             BoundaryValueImp,
                                             RangeImp >           derived_type;
}; // class AdvectionGodunovWithReconstructionLocalizableTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionGodunovWithReconstructionTraits
    : public AdvectionLaxFriedrichsTraits< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueImp, FVSpaceImp >
{
public:
  typedef AdvectionGodunovWithReconstruction< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueImp, FVSpaceImp > derived_type;
}; // class AdvectionLaxWendroffTraits

} // namespace internal


template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionLaxFriedrichsLocalizable
  : public Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp,
                                                                                LocalizableFunctionImp,
                                                                                SourceImp,
                                                                                BoundaryValueImp,
                                                                                RangeImp > >
  , public SystemAssembler< typename RangeImp::SpaceType >
{
  typedef Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp,
                                                                                LocalizableFunctionImp,
                                                                                SourceImp,
                                                                                BoundaryValueImp,
                                                                                RangeImp > >    OperatorBaseType;
  typedef SystemAssembler< typename RangeImp::SpaceType >                                       AssemblerBaseType;
public:
  typedef internal::AdvectionLaxFriedrichsLocalizableTraits< AnalyticalFluxImp,
                                                             LocalizableFunctionImp,
                                                             SourceImp,
                                                             BoundaryValueImp,
                                                             RangeImp >                       Traits;

  typedef typename Traits::GridViewType                                                       GridViewType;
  typedef typename Traits::SourceType                                                         SourceType;
  typedef typename Traits::RangeType                                                          RangeType;
  typedef typename Traits::AnalyticalFluxType                                                 AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType                                            LocalizableFunctionType;
  typedef typename Traits::BoundaryValueType::ExpressionFunctionType                          BoundaryValueType;

  typedef typename Dune::GDT::LocalEvaluation::LaxFriedrichs::Inner< LocalizableFunctionImp > NumericalFluxType;
  typedef typename Dune::GDT::LocalEvaluation::LaxFriedrichs::Dirichlet< LocalizableFunctionImp,
                                                                         BoundaryValueType >  NumericalBoundaryFluxType;
  typedef typename Dune::GDT::LocalOperator::Codim1FV< NumericalFluxType >                    LocalOperatorType;
  typedef typename Dune::GDT::LocalOperator::Codim1FVBoundary< NumericalBoundaryFluxType >    LocalBoundaryOperatorType;
  typedef typename LocalAssembler::Codim1CouplingFV< LocalOperatorType >                      InnerAssemblerType;
  typedef typename LocalAssembler::Codim1BoundaryFV< LocalBoundaryOperatorType >              BoundaryAssemblerType;

  AdvectionLaxFriedrichsLocalizable(const AnalyticalFluxType& analytical_flux,
                                    const LocalizableFunctionType& ratio_dt_dx,
                                    const SourceType& source,
                                    const BoundaryValueType& boundary_values,
                                    RangeType& range,
                                    const bool use_local)
    : OperatorBaseType()
    , AssemblerBaseType(range.space())
    , local_operator_(analytical_flux, ratio_dt_dx, use_local)
    , local_boundary_operator_(analytical_flux, ratio_dt_dx, boundary_values, use_local)
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



template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionLaxFriedrichs
  : public Dune::GDT::OperatorInterface< internal::AdvectionLaxFriedrichsTraits<  AnalyticalFluxImp,
                                                                                  LocalizableFunctionImp,
                                                                                  BoundaryValueImp,
                                                                                  FVSpaceImp > >
{
  typedef Dune::GDT::OperatorInterface< internal::AdvectionLaxFriedrichsTraits<  AnalyticalFluxImp,
                                                                                 LocalizableFunctionImp,
                                                                                 BoundaryValueImp,
                                                                                 FVSpaceImp > > OperatorBaseType;

public:
  typedef internal::AdvectionLaxFriedrichsTraits< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueImp, FVSpaceImp > Traits;
  typedef typename Traits::GridViewType            GridViewType;
  typedef typename Traits::AnalyticalFluxType      AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::BoundaryValueType       BoundaryValueType;
  typedef typename Traits::FVSpaceType             FVSpaceType;

  AdvectionLaxFriedrichs(const AnalyticalFluxType& analytical_flux,
                         const LocalizableFunctionType& ratio_dt_dx,
                         const BoundaryValueType& boundary_values,
                         const FVSpaceType& fv_space,
                         const bool use_local = false)
    : OperatorBaseType()
    , analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
    , boundary_values_(boundary_values)
    , fv_space_(fv_space)
    , use_local_(use_local)
  {}

  const GridViewType& grid_view() const
  {
    return fv_space_.grid_view();
  }

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range, const double time = 0.0) const
  {
    typename BoundaryValueType::ExpressionFunctionType current_boundary_values = boundary_values_.evaluate_at_time(time);
    AdvectionLaxFriedrichsLocalizable< AnalyticalFluxType,
                                       LocalizableFunctionType,
                                       SourceType,
                                       BoundaryValueType,
                                       RangeType > localizable_operator(analytical_flux_, ratio_dt_dx_, source, current_boundary_values, range, use_local_);
    localizable_operator.apply();
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
  const BoundaryValueType&  boundary_values_;
  const FVSpaceType& fv_space_;
  const bool use_local_;
}; // class AdvectionLaxFriedrichs

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionGodunovLocalizable
  : public Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionGodunovLocalizableTraits< AnalyticalFluxImp,
                                                                          LocalizableFunctionImp,
                                                                          SourceImp,
                                                                          BoundaryValueImp,
                                                                          RangeImp > >
  , public SystemAssembler< typename RangeImp::SpaceType >
{
  typedef Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionGodunovLocalizableTraits< AnalyticalFluxImp,
                                                                          LocalizableFunctionImp,
                                                                          SourceImp,
                                                                          BoundaryValueImp,
                                                                          RangeImp > >        OperatorBaseType;
  typedef SystemAssembler< typename RangeImp::SpaceType >                                     AssemblerBaseType;
public:
  typedef internal::AdvectionGodunovLocalizableTraits< AnalyticalFluxImp,
                                                       LocalizableFunctionImp,
                                                       SourceImp,
                                                       BoundaryValueImp,
                                                       RangeImp >                             Traits;

  typedef typename Traits::GridViewType                                                       GridViewType;
  typedef typename Traits::SourceType                                                         SourceType;
  typedef typename Traits::RangeType                                                          RangeType;
  typedef typename Traits::AnalyticalFluxType                                                 AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType                                            LocalizableFunctionType;
  typedef typename Traits::BoundaryValueType                                                  BoundaryValueType;

  typedef typename Dune::GDT::LocalEvaluation::Godunov::Inner< LocalizableFunctionImp >       NumericalFluxType;
  typedef typename Dune::GDT::LocalEvaluation::Godunov::Dirichlet< LocalizableFunctionImp,
                                                                   BoundaryValueType >        NumericalBoundaryFluxType;
  typedef typename Dune::GDT::LocalOperator::Codim1FV< NumericalFluxType >                    LocalOperatorType;
  typedef typename Dune::GDT::LocalOperator::Codim1FVBoundary< NumericalBoundaryFluxType >    LocalBoundaryOperatorType;
  typedef typename LocalAssembler::Codim1CouplingFV< LocalOperatorType >                      InnerAssemblerType;
  typedef typename LocalAssembler::Codim1BoundaryFV< LocalBoundaryOperatorType >              BoundaryAssemblerType;

  AdvectionGodunovLocalizable(const AnalyticalFluxType& analytical_flux,
                              const LocalizableFunctionType& ratio_dt_dx,
                              const SourceType& source,
                              const BoundaryValueType& boundary_values,
                              RangeType& range,
                              const bool is_linear)
    : OperatorBaseType()
    , AssemblerBaseType(range.space())
    , local_operator_(analytical_flux, ratio_dt_dx, is_linear)
    , local_boundary_operator_(analytical_flux, ratio_dt_dx, boundary_values, is_linear)
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
}; // class AdvectionGodunovLocalizable



template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionGodunov
  : public Dune::GDT::OperatorInterface< internal::AdvectionGodunovTraits<  AnalyticalFluxImp,
                                                                            LocalizableFunctionImp,
                                                                            BoundaryValueImp,
                                                                            FVSpaceImp > >
{
  typedef Dune::GDT::OperatorInterface< internal::AdvectionGodunovTraits<  AnalyticalFluxImp,
                                                                           LocalizableFunctionImp,
                                                                           BoundaryValueImp,
                                                                           FVSpaceImp > > OperatorBaseType;

public:
  typedef internal::AdvectionGodunovTraits< AnalyticalFluxImp,
                                            LocalizableFunctionImp,
                                            BoundaryValueImp,
                                            FVSpaceImp >          Traits;
  typedef typename Traits::GridViewType                           GridViewType;
  typedef typename Traits::AnalyticalFluxType                     AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType                LocalizableFunctionType;
  typedef typename Traits::BoundaryValueType                      BoundaryValueType;
  typedef typename Traits::FVSpaceType                            FVSpaceType;

  AdvectionGodunov(const AnalyticalFluxType& analytical_flux,
                         const LocalizableFunctionType& ratio_dt_dx,
                         const BoundaryValueType& boundary_values,
                         const FVSpaceType& fv_space,
                         const bool is_linear = false)
    : OperatorBaseType()
    , analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
    , boundary_values_(boundary_values)
    , fv_space_(fv_space)
    , is_linear_(is_linear)
  {}

  const GridViewType& grid_view() const
  {
    return fv_space_.grid_view();
  }

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range, const double time = 0.0) const
  {
    typename BoundaryValueType::ExpressionFunctionType current_boundary_values = boundary_values_.evaluate_at_time(time);
    AdvectionGodunovLocalizable<       AnalyticalFluxType,
                                       LocalizableFunctionType,
                                       SourceType,
                                       typename BoundaryValueType::ExpressionFunctionType,
                                       RangeType > localizable_operator(analytical_flux_, ratio_dt_dx_, source, current_boundary_values, range, is_linear_);
    localizable_operator.apply();
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
  const BoundaryValueType&  boundary_values_;
  const FVSpaceType& fv_space_;
  const bool is_linear_;
}; // class AdvectionGodunov

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionLaxWendroffLocalizable
  : public Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionLaxWendroffLocalizableTraits< AnalyticalFluxImp,
                                                                          LocalizableFunctionImp,
                                                                          SourceImp,
                                                                          BoundaryValueImp,
                                                                          RangeImp > >
  , public SystemAssembler< typename RangeImp::SpaceType >
{
  typedef Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionLaxWendroffLocalizableTraits< AnalyticalFluxImp,
                                                                          LocalizableFunctionImp,
                                                                          SourceImp,
                                                                          BoundaryValueImp,
                                                                          RangeImp > >        OperatorBaseType;
  typedef SystemAssembler< typename RangeImp::SpaceType >                                     AssemblerBaseType;
public:
  typedef internal::AdvectionLaxWendroffLocalizableTraits< AnalyticalFluxImp,
                                                       LocalizableFunctionImp,
                                                       SourceImp,
                                                       BoundaryValueImp,
                                                       RangeImp >                             Traits;

  typedef typename Traits::GridViewType                                                       GridViewType;
  typedef typename Traits::SourceType                                                         SourceType;
  typedef typename Traits::RangeType                                                          RangeType;
  typedef typename Traits::AnalyticalFluxType                                                 AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType                                            LocalizableFunctionType;
  typedef typename Traits::BoundaryValueType                                                  BoundaryValueType;

  typedef typename Dune::GDT::LocalEvaluation::LaxWendroff::Inner< LocalizableFunctionImp >       NumericalFluxType;
  typedef typename Dune::GDT::LocalEvaluation::LaxWendroff::Dirichlet< LocalizableFunctionImp,
                                                                   BoundaryValueType >        NumericalBoundaryFluxType;
  typedef typename Dune::GDT::LocalOperator::Codim1FV< NumericalFluxType >                    LocalOperatorType;
  typedef typename Dune::GDT::LocalOperator::Codim1FVBoundary< NumericalBoundaryFluxType >    LocalBoundaryOperatorType;
  typedef typename LocalAssembler::Codim1CouplingFV< LocalOperatorType >                      InnerAssemblerType;
  typedef typename LocalAssembler::Codim1BoundaryFV< LocalBoundaryOperatorType >              BoundaryAssemblerType;

  AdvectionLaxWendroffLocalizable(const AnalyticalFluxType& analytical_flux,
                              const LocalizableFunctionType& ratio_dt_dx,
                              const SourceType& source,
                              const BoundaryValueType& boundary_values,
                              RangeType& range,
                              const bool is_linear)
    : OperatorBaseType()
    , AssemblerBaseType(range.space())
    , local_operator_(analytical_flux, ratio_dt_dx, is_linear)
    , local_boundary_operator_(analytical_flux, ratio_dt_dx, boundary_values, is_linear)
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
}; // class AdvectionLaxWendroffLocalizable



template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionLaxWendroff
  : public Dune::GDT::OperatorInterface< internal::AdvectionLaxWendroffTraits<  AnalyticalFluxImp,
                                                                            LocalizableFunctionImp,
                                                                            BoundaryValueImp,
                                                                            FVSpaceImp > >
{
  typedef Dune::GDT::OperatorInterface< internal::AdvectionLaxWendroffTraits<  AnalyticalFluxImp,
                                                                           LocalizableFunctionImp,
                                                                           BoundaryValueImp,
                                                                           FVSpaceImp > > OperatorBaseType;

public:
  typedef internal::AdvectionLaxWendroffTraits< AnalyticalFluxImp,
                                            LocalizableFunctionImp,
                                            BoundaryValueImp,
                                            FVSpaceImp >          Traits;
  typedef typename Traits::GridViewType                           GridViewType;
  typedef typename Traits::AnalyticalFluxType                     AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType                LocalizableFunctionType;
  typedef typename Traits::BoundaryValueType                      BoundaryValueType;
  typedef typename Traits::FVSpaceType                            FVSpaceType;

  AdvectionLaxWendroff(const AnalyticalFluxType& analytical_flux,
                         const LocalizableFunctionType& ratio_dt_dx,
                         const BoundaryValueType& boundary_values,
                         const FVSpaceType& fv_space,
                         const bool is_linear = false)
    : OperatorBaseType()
    , analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
    , boundary_values_(boundary_values)
    , fv_space_(fv_space)
    , is_linear_(is_linear)
  {}

  const GridViewType& grid_view() const
  {
    return fv_space_.grid_view();
  }

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range, const double time = 0.0) const
  {
    typename BoundaryValueType::ExpressionFunctionType current_boundary_values = boundary_values_.evaluate_at_time(time);
    AdvectionLaxWendroffLocalizable<   AnalyticalFluxType,
                                       LocalizableFunctionType,
                                       SourceType,
                                       typename BoundaryValueType::ExpressionFunctionType,
                                       RangeType > localizable_operator(analytical_flux_, ratio_dt_dx_, source, current_boundary_values, range, is_linear_);
    localizable_operator.apply();
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
  const BoundaryValueType&  boundary_values_;
  const FVSpaceType& fv_space_;
  const bool is_linear_;
}; // class AdvectionLaxWendroff

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
class AdvectionGodunovWithReconstructionLocalizable
  : public Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionGodunovWithReconstructionLocalizableTraits< AnalyticalFluxImp,
                                                                          LocalizableFunctionImp,
                                                                          SourceImp,
                                                                          BoundaryValueImp,
                                                                          RangeImp > >
  , public SystemAssembler< typename RangeImp::SpaceType >
{
  typedef Dune::GDT::LocalizableOperatorInterface<
                             internal::AdvectionGodunovWithReconstructionLocalizableTraits< AnalyticalFluxImp,
                                                                          LocalizableFunctionImp,
                                                                          SourceImp,
                                                                          BoundaryValueImp,
                                                                          RangeImp > >        OperatorBaseType;
  typedef SystemAssembler< typename RangeImp::SpaceType >                                     AssemblerBaseType;
public:
  typedef internal::AdvectionGodunovWithReconstructionLocalizableTraits< AnalyticalFluxImp,
                                                       LocalizableFunctionImp,
                                                       SourceImp,
                                                       BoundaryValueImp,
                                                       RangeImp >                             Traits;

  typedef typename Traits::GridViewType                                                       GridViewType;
  typedef typename Traits::SourceType                                                         SourceType;
  typedef typename Traits::RangeType                                                          RangeType;
  typedef typename Traits::RangeFieldType                                                     RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                                                 AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType                                            LocalizableFunctionType;
  typedef typename Traits::BoundaryValueType                                                  BoundaryValueType;

  typedef typename Dune::GDT::LocalEvaluation::Godunov::Inner< LocalizableFunctionImp >       NumericalFluxType;
  typedef typename Dune::GDT::LocalEvaluation::Godunov::Dirichlet< LocalizableFunctionImp,
                                                                   BoundaryValueType >        NumericalBoundaryFluxType;
  typedef typename Dune::GDT::LocalOperator::Codim1FV< NumericalFluxType >                    LocalOperatorType;
  typedef typename Dune::GDT::LocalOperator::Codim1FVBoundary< NumericalBoundaryFluxType >    LocalBoundaryOperatorType;
  typedef typename LocalAssembler::Codim1CouplingFV< LocalOperatorType >                      InnerAssemblerType;
  typedef typename LocalAssembler::Codim1BoundaryFV< LocalBoundaryOperatorType >              BoundaryAssemblerType;

  typedef typename GridViewType::Grid                                       GridType;
  typedef typename GridType::template Codim< 0 >::EntityPointer             EntityPointerType;
  static const size_t dimRange = SourceType::dimRange;
  static const size_t dimRangeCols = SourceType::dimRangeCols;

  typedef Spaces::DG::PdelabBasedProduct< GridViewType,
                                          1, //polOrder
                                          RangeFieldType,
                                          dimRange,
                                          dimRangeCols >               DGSpaceType;
  typedef DiscreteFunction< DGSpaceType, typename SourceType::VectorType >  ReconstructedDiscreteFunctionType;

  AdvectionGodunovWithReconstructionLocalizable(const AnalyticalFluxType& analytical_flux,
                                                const LocalizableFunctionType& dx,
                                                const double dt,
                                                const SourceType& source,
                                                const BoundaryValueType& boundary_values,
                                                RangeType& range,
                                                const bool is_linear)
    : OperatorBaseType()
    , AssemblerBaseType(range.space())
    , dx_(dx)
    , local_operator_(analytical_flux, dx, dt, is_linear)
    , boundary_values_(boundary_values)
    , local_boundary_operator_(analytical_flux, dx, dt, boundary_values_, is_linear)
    , inner_assembler_(local_operator_)
    , boundary_assembler_(local_boundary_operator_)
    , source_(source)
    , range_(range)
    , grid_view_(source.space().grid_view())
    , dg_space_(grid_view_)
    , reconstruction_(dg_space_, "reconstructed")
  {
    reconstruct_linear();
  }

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
    const ReconstructedDiscreteFunctionType reconstruction(reconstruction_);
    this->add(inner_assembler_, reconstruction, range_, new DSG::ApplyOn::InnerIntersections< GridViewType >());
    this->add(inner_assembler_, reconstruction, range_, new DSG::ApplyOn::PeriodicIntersections< GridViewType >());
    this->add(boundary_assembler_, reconstruction, range_, new DSG::ApplyOn::NonPeriodicBoundaryIntersections< GridViewType >());
    this->assemble();
  }

private:
  void reconstruct_linear() {
    // take slope sigma_i, then reconstruct linear on entity x as u_i + sigma_i*(x-(x+dx/2))
    // possible slopes:

    // walk the grid to reconstruct
    const auto it_end = grid_view_.template end< 0 >();
    for (auto it = grid_view_.template begin< 0 >(); it != it_end; ++it) {
      const auto& entity = *it;
      // create EntityPointers for current entity and two neighbors to the right and left
      EntityPointerType entity_ptr(entity);
      EntityPointerType left_neighbor_ptr(entity);
      EntityPointerType second_left_neighbor_ptr(entity);
      EntityPointerType right_neighbor_ptr(entity);
      EntityPointerType second_right_neighbor_ptr(entity);
      typename SourceType::RangeType right_boundary_value;
      typename SourceType::RangeType left_boundary_value;
      // get local_discrete_function of reconstruction
      bool on_left_boundary(false);
      bool on_right_boundary(false);
      bool left_neighbor_on_boundary(false);
      bool right_neighbor_on_boundary(false);
      // walk over intersections to get neighbors
      const auto i_it_end = grid_view_.iend(entity);
      for (auto i_it = grid_view_.ibegin(entity); i_it != i_it_end; ++i_it) {
        const auto& intersection = *i_it;
        const auto& neighbor = intersection.neighbor() ? *(intersection.outside()) : entity;
        if (intersection.neighbor()) {
          if ((neighbor.geometry().center()[0] < entity.geometry().center()[0] && !(intersection.boundary())) || (neighbor.geometry().center()[0] > entity.geometry().center()[0] && intersection.boundary())) {
            left_neighbor_ptr = EntityPointerType(neighbor);
            const auto left_neighbor_i_it_end = grid_view_.iend(neighbor);
            for (auto left_neighbor_i_it = grid_view_.ibegin(neighbor); left_neighbor_i_it != left_neighbor_i_it_end; ++left_neighbor_i_it) {
              const auto& neighbor_intersection = *left_neighbor_i_it;
              const auto& neighbor_neighbor = neighbor_intersection.neighbor() ? *(neighbor_intersection.outside()) : entity;
              if (neighbor_intersection.neighbor()) {
                if ((neighbor_neighbor.geometry().center()[0] < neighbor.geometry().center()[0] && !(neighbor_intersection.boundary())) || (neighbor_neighbor.geometry().center()[0] > neighbor.geometry().center()[0] && neighbor_intersection.boundary()))
                  second_left_neighbor_ptr = EntityPointerType(neighbor_neighbor);
              } else {
                if (neighbor_intersection.geometry().center()[0] < neighbor.geometry().center()[0]) {
                  left_neighbor_on_boundary = true;
                  left_boundary_value = boundary_values_.local_function(neighbor)->evaluate(neighbor.geometry().local(neighbor_intersection.geometry().center()));
                } else {
                  DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
                }
              }
            }
          } else {
            right_neighbor_ptr = EntityPointerType(neighbor);
            const auto right_neighbor_i_it_end = grid_view_.iend(neighbor);
            for (auto right_neighbor_i_it = grid_view_.ibegin(neighbor); right_neighbor_i_it != right_neighbor_i_it_end; ++right_neighbor_i_it) {
              const auto& neighbor_intersection = *right_neighbor_i_it;
              const auto& neighbor_neighbor = neighbor_intersection.neighbor() ? *(neighbor_intersection.outside()) : entity;
              if (neighbor_intersection.neighbor()) {
                if ((neighbor_neighbor.geometry().center()[0] > neighbor.geometry().center()[0] && !(neighbor_intersection.boundary())) || (neighbor_neighbor.geometry().center()[0] < neighbor.geometry().center()[0] && neighbor_intersection.boundary()))
                  second_right_neighbor_ptr = EntityPointerType(neighbor_neighbor);
              } else {
                if (neighbor_intersection.geometry().center()[0] < neighbor.geometry().center()[0]) {
                  DUNE_THROW(Dune::InvalidStateException, "This should not happen!");
                } else {
                  right_neighbor_on_boundary = true;
                  right_boundary_value = boundary_values_.local_function(neighbor)->evaluate(neighbor.geometry().local(neighbor_intersection.geometry().center()));
                }
              }
            }
          }
        } else {
          if (intersection.geometry().center()[0] < entity.geometry().center()[0]) {
            on_left_boundary = true;
            left_boundary_value = boundary_values_.local_function(entity)->evaluate(entity.geometry().local(intersection.geometry().center()));
          } else {
            on_right_boundary = true;
            right_boundary_value = boundary_values_.local_function(entity)->evaluate(entity.geometry().local(intersection.geometry().center()));
          }
        }
      }
      // get values of discrete function on the five point stencil
      const auto u_second_left = (on_left_boundary || left_neighbor_on_boundary) ? left_boundary_value : source_.local_discrete_function(*second_left_neighbor_ptr)->evaluate(second_left_neighbor_ptr->geometry().local(second_left_neighbor_ptr->geometry().center()));
      const auto u_left = on_left_boundary ? left_boundary_value : source_.local_discrete_function(*left_neighbor_ptr)->evaluate(left_neighbor_ptr->geometry().local(left_neighbor_ptr->geometry().center()));
      const auto u_second_right = (on_right_boundary || right_neighbor_on_boundary) ? right_boundary_value : source_.local_discrete_function(*second_right_neighbor_ptr)->evaluate(second_right_neighbor_ptr->geometry().local(second_right_neighbor_ptr->geometry().center()));
      const auto u_right = on_right_boundary ? right_boundary_value : source_.local_discrete_function(*right_neighbor_ptr)->evaluate(right_neighbor_ptr->geometry().local(right_neighbor_ptr->geometry().center()));
      const auto u_entity = source_.local_discrete_function(*entity_ptr)->evaluate(entity_ptr->geometry().local(entity_ptr->geometry().center()));
      // calculate slope
      const double dx = (dx_.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center())))[0];
      for (size_t factor_index = 0; factor_index < dimRange; ++factor_index) {
        RangeFieldType slope_left = u_entity[factor_index] - u_left[factor_index];
        slope_left *= 1.0/dx;
        RangeFieldType slope_right = u_right[factor_index] - u_entity[factor_index];
        slope_right *= 1.0/dx;
        RangeFieldType centered_slope = u_right[factor_index] - u_left[factor_index];
        centered_slope *= 1.0/(2.0*dx);
        RangeFieldType slope = minmod(slope_left, slope_right, centered_slope);
        RangeFieldType local_slope = slope*dx;
        //      std::cout << "entity: " << DSC::toString(entity.geometry().center()) << " and slope: " << local_slope << std::endl;
        //      std::cout << "u_left: " << u_left << " and u_right" << u_right << std::endl;
        //      std::cout << "left neighbor: " << DSC::toString(left_neighbor_ptr->geometry().center()) << " and right " << DSC::toString(right_neighbor_ptr->geometry().center()) << std::endl;
        // calculate value of the reconstruction on left and right side of the (reference) entity
        const RangeFieldType value_left = u_entity[factor_index] - 0.5*local_slope;
        const RangeFieldType value_right = u_entity[factor_index] + 0.5*local_slope;
        // set values on dofs, dof with local index 0 for each factor space corresponds to 1 - x, local index 1 to x
        reconstruction_.vector().set_entry(reconstruction_.space().factor_mapper().mapToGlobal(factor_index, entity, 0), value_left);
        reconstruction_.vector().set_entry(reconstruction_.space().factor_mapper().mapToGlobal(factor_index, entity, 1), value_right);
//        if (factor_index == 0)
//          reconstruction_.template visualize_factor< 0 >("reconstructed" + DSC::toString(step_number), false);
      }
    }
    ++step_number;
  }

  RangeFieldType minmod(const RangeFieldType slope_left, const RangeFieldType slope_right, const RangeFieldType /*centered_slope*/ = RangeFieldType(0)) const
  {
    RangeFieldType slope_left_abs = std::abs(slope_left);
    RangeFieldType slope_right_abs = std::abs(slope_right);
    if (slope_left_abs < slope_right_abs && slope_left*slope_right > 0)
      return slope_left;
    else if (DSC::FloatCmp::ge(slope_left_abs, slope_right_abs) && slope_left*slope_right > 0)
      return slope_right;
    else
      return 0.0;
  }

  RangeFieldType maxmod(const RangeFieldType slope_left, const RangeFieldType slope_right, const RangeFieldType /*centered_slope*/ = RangeFieldType(0)) const
  {
    RangeFieldType slope_left_abs = std::abs(slope_left);
    RangeFieldType slope_right_abs = std::abs(slope_right);
    if (slope_left_abs > slope_right_abs && slope_left*slope_right > 0)
      return slope_left;
    else if (DSC::FloatCmp::le(slope_left_abs, slope_right_abs) && slope_left*slope_right > 0)
      return slope_right;
    else
      return 0.0;
  }

  RangeFieldType superbee(const RangeFieldType slope_left, const RangeFieldType slope_right, const RangeFieldType /*centered_slope*/ = RangeFieldType(0)) const
  {
    return maxmod(minmod(slope_left, 2.0*slope_right), minmod(2.0*slope_left, slope_right));
  }

  RangeFieldType mc(const RangeFieldType slope_left, const RangeFieldType slope_right, const RangeFieldType centered_slope) const
  {
    return minmod(minmod(2.0*slope_left, 2.0*slope_right), centered_slope);
  }

  const LocalizableFunctionType& dx_;
  const LocalOperatorType local_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
  const InnerAssemblerType inner_assembler_;
  const BoundaryValueType boundary_values_;
  const BoundaryAssemblerType boundary_assembler_;
  const SourceType& source_;
  RangeType& range_;
  const GridViewType& grid_view_;
  const DGSpaceType dg_space_;
  ReconstructedDiscreteFunctionType reconstruction_;
  static size_t step_number;
}; // class AdvectionGodunovWithReconstructionLocalizable

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class SourceImp, class BoundaryValueImp, class RangeImp >
size_t AdvectionGodunovWithReconstructionLocalizable< AnalyticalFluxImp, LocalizableFunctionImp, SourceImp, BoundaryValueImp, RangeImp >::step_number = 0;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueImp, class FVSpaceImp >
class AdvectionGodunovWithReconstruction
  : public Dune::GDT::OperatorInterface< internal::AdvectionGodunovWithReconstructionTraits<  AnalyticalFluxImp,
                                                                            LocalizableFunctionImp,
                                                                            BoundaryValueImp,
                                                                            FVSpaceImp > >
{
  typedef Dune::GDT::OperatorInterface< internal::AdvectionGodunovWithReconstructionTraits<  AnalyticalFluxImp,
                                                                           LocalizableFunctionImp,
                                                                           BoundaryValueImp,
                                                                           FVSpaceImp > > OperatorBaseType;

public:
  typedef internal::AdvectionGodunovWithReconstructionTraits< AnalyticalFluxImp,
                                            LocalizableFunctionImp,
                                            BoundaryValueImp,
                                            FVSpaceImp >          Traits;
  typedef typename Traits::GridViewType                           GridViewType;
  typedef typename Traits::AnalyticalFluxType                     AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType                LocalizableFunctionType;
  typedef typename Traits::BoundaryValueType                      BoundaryValueType;
  typedef typename Traits::FVSpaceType                            FVSpaceType;

  AdvectionGodunovWithReconstruction(const AnalyticalFluxType& analytical_flux,
                         const LocalizableFunctionType& dx,
                         const double dt,
                         const BoundaryValueType& boundary_values,
                         const FVSpaceType& fv_space,
                         const bool is_linear = false)
    : OperatorBaseType()
    , analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , boundary_values_(boundary_values)
    , fv_space_(fv_space)
    , is_linear_(is_linear)
  {}

  const GridViewType& grid_view() const
  {
    return fv_space_.grid_view();
  }

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range, const double time = 0.0) const
  {
    typename BoundaryValueType::ExpressionFunctionType current_boundary_values = boundary_values_.evaluate_at_time(time);
    AdvectionGodunovWithReconstructionLocalizable< AnalyticalFluxType,
                                       LocalizableFunctionType,
                                       SourceType,
                                       typename BoundaryValueType::ExpressionFunctionType,
                                       RangeType > localizable_operator(analytical_flux_, dx_, dt_, source, current_boundary_values, range, is_linear_);
    localizable_operator.apply();
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const BoundaryValueType& boundary_values_;
  const FVSpaceType& fv_space_;
  const bool is_linear_;
}; // class AdvectionGodunovWithReconstruction

} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_HH
