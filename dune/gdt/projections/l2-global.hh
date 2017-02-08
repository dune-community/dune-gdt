// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_GDT_PROJECTIONS_L2_GLOBAL_HH
#define DUNE_GDT_PROJECTIONS_L2_GLOBAL_HH

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/solver.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/operators/l2.hh>

namespace Dune {
namespace GDT {


// forward
template <class GridViewImp, class FieldImp = double>
class L2GlobalProjectionOperator;


namespace internal {


template <class GridViewImp, class FieldImp>
class L2GlobalProjectionOperatorTraits
{
public:
  typedef L2GlobalProjectionOperator<GridViewImp, FieldImp> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class GridViewImp, class SourceImp, class RangeImp>
class L2GlobalProjectionLocalizableOperator : public LocalizableOperatorBase<GridViewImp, SourceImp, RangeImp>
{
  static_assert(XT::Functions::is_localizable_function<SourceImp>::value, "");
  static_assert(is_discrete_function<RangeImp>::value, "");
  typedef LocalizableOperatorBase<GridViewImp, SourceImp, RangeImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::SourceType;
  using typename BaseType::RangeType;
  typedef typename RangeType::VectorType VectorType;
  typedef
      typename XT::LA::Container<typename RangeType::RangeFieldType, VectorType::Traits::sparse_matrix_type>::MatrixType

          MatrixType;

private:
  typedef typename RangeType::SpaceType SpaceType;
  typedef L2MatrixOperator<SpaceType, MatrixType, GridViewType> LhsOperatorType;
  typedef L2VolumeVectorFunctional<SourceType, SpaceType, VectorType, GridViewType> RhsFunctionalType;

public:
  L2GlobalProjectionLocalizableOperator(const size_t over_integrate,
                                        GridViewType grd_vw,
                                        const SourceType& src,
                                        RangeType& rng)
    : BaseType(grd_vw, src, rng)
    , lhs_operator_(over_integrate, range_.space(), BaseType::grid_view())
    , rhs_functional_(over_integrate, source_, range_.space(), BaseType::grid_view())
    , solved_(false)
  {
    this->append(lhs_operator_);
    this->append(rhs_functional_);
    issue_warning(this->range().space());
  }

  L2GlobalProjectionLocalizableOperator(GridViewType grd_vw, const SourceType& src, RangeType& rng)
    : L2GlobalProjectionLocalizableOperator(0, grd_vw, src, rng)
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


template <class GridViewType, class SourceType, class SpaceType, class VectorType>
typename std::
    enable_if<XT::Grid::is_layer<GridViewType>::value && XT::Functions::is_localizable_function<SourceType>::value
                  && is_space<SpaceType>::value
                  && XT::LA::is_vector<VectorType>::value,
              std::unique_ptr<L2GlobalProjectionLocalizableOperator<GridViewType,
                                                                    SourceType,
                                                                    DiscreteFunction<SpaceType, VectorType>>>>::type
    make_global_l2_projection_localizable_operator(const GridViewType& grid_view,
                                                   const SourceType& source,
                                                   DiscreteFunction<SpaceType, VectorType>& range,
                                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2GlobalProjectionLocalizableOperator<GridViewType,
                                                                             SourceType,
                                                                             DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, grid_view, source, range);
} // ... make_global_l2_projection_localizable_operator(...)

template <class SourceType, class SpaceType, class VectorType>
typename std::
    enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                  && XT::LA::is_vector<VectorType>::value,
              std::unique_ptr<L2GlobalProjectionLocalizableOperator<typename SpaceType::GridViewType,
                                                                    SourceType,
                                                                    DiscreteFunction<SpaceType, VectorType>>>>::type
    make_global_l2_projection_localizable_operator(const SourceType& source,
                                                   DiscreteFunction<SpaceType, VectorType>& range,
                                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2GlobalProjectionLocalizableOperator<typename SpaceType::GridViewType,
                                                                             SourceType,
                                                                             DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, range.space().grid_view(), source, range);
} // ... make_global_l2_projection_localizable_operator(...)


template <class GridViewImp, class FieldImp>
class L2GlobalProjectionOperator
    : public OperatorInterface<internal::L2GlobalProjectionOperatorTraits<GridViewImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2GlobalProjectionOperatorTraits<GridViewImp, FieldImp>> BaseType;

public:
  typedef internal::L2GlobalProjectionOperatorTraits<GridViewImp, FieldImp> Traits;
  typedef GridViewImp GridViewType;
  using typename BaseType::FieldType;

private:
  typedef typename XT::Grid::Entity<GridViewType>::Type E;
  typedef typename GridViewType::ctype D;
  static const size_t d = GridViewType::dimension;

public:
  L2GlobalProjectionOperator(const size_t over_integrate, GridViewType grid_view)
    : grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {
  }

  L2GlobalProjectionOperator(GridViewType grid_view)
    : grid_view_(grid_view)
    , over_integrate_(0)
  {
  }

  template <class R, size_t r, size_t rC, class S, class V>
  void apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
             DiscreteFunction<S, V>& range) const
  {
    typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC> SourceType;
    L2GlobalProjectionLocalizableOperator<GridViewType, SourceType, DiscreteFunction<S, V>> op(
        over_integrate_, grid_view_, source, range);
    op.apply();
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/, const SourceType& /*source*/) const
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
  GridViewType grid_view_;
  const size_t over_integrate_;
}; // class L2GlobalProjectionOperator


template <class GridViewType>
typename std::enable_if<XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<L2GlobalProjectionOperator<GridViewType>>>::type
make_global_l2_projection_operator(const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2GlobalProjectionOperator<GridViewType>>(over_integrate, grid_view);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_L2_GLOBAL_HH
