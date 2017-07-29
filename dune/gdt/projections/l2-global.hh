// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_PROJECTIONS_L2_GLOBAL_HH
#define DUNE_GDT_PROJECTIONS_L2_GLOBAL_HH

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/operators/l2.hh>

namespace Dune {
namespace GDT {


// forward
template <class GridLayerImp, class FieldImp = double>
class L2GlobalProjectionOperator;


namespace internal {


template <class GridLayerImp, class FieldImp>
class L2GlobalProjectionOperatorTraits
{
public:
  typedef L2GlobalProjectionOperator<GridLayerImp, FieldImp> derived_type;
  typedef NoJacobian JacobianType;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class GridLayerImp, class SourceImp, class RangeImp>
class L2GlobalProjectionLocalizableOperator : public LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp>
{
  static_assert(XT::Functions::is_localizable_function<SourceImp>::value, "");
  static_assert(is_discrete_function<RangeImp>::value, "");
  typedef LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::SourceType;
  using typename BaseType::RangeType;
  typedef typename RangeType::VectorType VectorType;
  typedef
      typename XT::LA::Container<typename RangeType::RangeFieldType, VectorType::Traits::sparse_matrix_type>::MatrixType

          MatrixType;

private:
  typedef typename RangeType::SpaceType SpaceType;
  typedef L2MatrixOperator<SpaceType, MatrixType, GridLayerType> LhsOperatorType;
  typedef L2VolumeVectorFunctional<SourceType, SpaceType, VectorType, GridLayerType> RhsFunctionalType;

public:
  L2GlobalProjectionLocalizableOperator(const size_t over_integrate,
                                        GridLayerType grd_vw,
                                        const SourceType& src,
                                        RangeType& rng,
                                        const XT::Common::Parameter& param = {})
    : BaseType(grd_vw, src, rng)
    , lhs_operator_(over_integrate, range_.space(), BaseType::grid_layer())
    , rhs_functional_(over_integrate, source_, range_.space(), BaseType::grid_layer())
    , solved_(false)
  {
    this->append(lhs_operator_);
    this->append(rhs_functional_);
    issue_warning(this->range().space());
  }

  L2GlobalProjectionLocalizableOperator(GridLayerType grd_vw,
                                        const SourceType& src,
                                        RangeType& rng,
                                        const XT::Common::Parameter& param = {})
    : L2GlobalProjectionLocalizableOperator(0, grd_vw, src, rng, param)
  {
  }

  void apply()
  {
    if (solved_)
      return;
    BaseType::apply();
    try {
      XT::LA::Solver<MatrixType>(lhs_operator_.matrix()).apply(rhs_functional_.vector(), range_.vector());
    } catch (XT::Common::Exceptions::linear_solver_failed& ee) {
      DUNE_THROW(projection_error,
                 "L2 projection failed because a global matrix could not be inverted!\n\n"
                     << "This was the original error: "
                     << ee.what());
    }
  } // ... apply(...)

private:
  template <class S, bool warn>
  struct Warning
  {
    static inline void issue()
    {
    }
  };

  template <class S>
  struct Warning<S, true>
  {
    static inline void issue()
    {
      Dune::XT::Common::TimedLogger().get("gdt.l2globalprojectionlocalizableoperator").warn()
          << "You are using this operator to project onto a discontinuous discrete function space (see below)!\n"
          << "Consider to use L2LocalProjectionLocalizableOperator instead!\n"
          << "You can disable this warning by defining "
          << "DUNE_GDT_PROJECTIONS_L2_GLOBAL_LOCALIZABLE_DISABLE_WARNING\n"
          << "at compile time or by disabling the Dune::XT::Common::TimedLogger() instance at runtime.\n"
          << "The type of the range space is: " << Dune::XT::Common::Typename<S>::value() << std::endl;
    } // ... issue_warning(...)
  };

  template <class S>
  static inline void issue_warning(const S&)
  {
#ifndef DUNE_GDT_PROJECTIONS_L2_GLOBAL_LOCALIZABLE_DISABLE_WARNING
    Warning<S, !S::continuous>::issue();
#endif
  }

  using BaseType::source_;
  using BaseType::range_;

  LhsOperatorType lhs_operator_;
  RhsFunctionalType rhs_functional_;
  bool solved_;
}; // class L2GlobalProjectionLocalizableOperator


template <class GridLayerType, class SourceType, class SpaceType, class VectorType>
typename std::
    enable_if<XT::Grid::is_layer<GridLayerType>::value && XT::Functions::is_localizable_function<SourceType>::value
                  && is_space<SpaceType>::value
                  && XT::LA::is_vector<VectorType>::value,
              std::unique_ptr<L2GlobalProjectionLocalizableOperator<GridLayerType,
                                                                    SourceType,
                                                                    DiscreteFunction<SpaceType, VectorType>>>>::type
    make_global_l2_projection_localizable_operator(const GridLayerType& grid_layer,
                                                   const SourceType& source,
                                                   DiscreteFunction<SpaceType, VectorType>& range,
                                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2GlobalProjectionLocalizableOperator<GridLayerType,
                                                                             SourceType,
                                                                             DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, grid_layer, source, range);
} // ... make_global_l2_projection_localizable_operator(...)

template <class SourceType, class SpaceType, class VectorType>
typename std::
    enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                  && XT::LA::is_vector<VectorType>::value,
              std::unique_ptr<L2GlobalProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                                                    SourceType,
                                                                    DiscreteFunction<SpaceType, VectorType>>>>::type
    make_global_l2_projection_localizable_operator(const SourceType& source,
                                                   DiscreteFunction<SpaceType, VectorType>& range,
                                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2GlobalProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                                                             SourceType,
                                                                             DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, range.space().grid_layer(), source, range);
} // ... make_global_l2_projection_localizable_operator(...)


template <class GridLayerImp, class FieldImp>
class L2GlobalProjectionOperator
    : public OperatorInterface<internal::L2GlobalProjectionOperatorTraits<GridLayerImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2GlobalProjectionOperatorTraits<GridLayerImp, FieldImp>> BaseType;

public:
  typedef internal::L2GlobalProjectionOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::FieldType;

private:
  using E = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype D;
  static const size_t d = GridLayerType::dimension;

public:
  L2GlobalProjectionOperator(const size_t over_integrate, GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {
  }

  L2GlobalProjectionOperator(GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(0)
  {
  }

  template <class R, size_t r, size_t rC, class S, class V>
  void apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
             DiscreteFunction<S, V>& range,
             const XT::Common::Parameter& param = {}) const
  {
    typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC> SourceType;
    L2GlobalProjectionLocalizableOperator<GridLayerType, SourceType, DiscreteFunction<S, V>> op(
        over_integrate_, grid_layer_, source, range, param);
    op.apply();
  }

  template <class RangeType, class SourceType>
  FieldType
  apply2(const RangeType& /*range*/, const SourceType& /*source*/, const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  template <class RangeType, class SourceType>
  void
  apply_inverse(const RangeType& /*range*/, SourceType& /*source*/, const XT::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

private:
  GridLayerType grid_layer_;
  const size_t over_integrate_;
}; // class L2GlobalProjectionOperator


template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<L2GlobalProjectionOperator<GridLayerType>>>::type
make_global_l2_projection_operator(const GridLayerType& grid_layer, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2GlobalProjectionOperator<GridLayerType>>(over_integrate, grid_layer);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_L2_GLOBAL_HH
