// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_L2_GLOBAL_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_L2_GLOBAL_HH

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/operators/default.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/operators/l2.hh>

#include "../default.hh"

namespace Dune {
namespace GDT {


template <class GridViewImp, class SourceImp, class RangeImp,
          Stuff::LA::ChooseBackend container_backend = Stuff::LA::default_backend>
class L2GlobalProjectionLocalizableOperator : public LocalizableOperatorDefault<GridViewImp, SourceImp, RangeImp>
{
  static_assert(Stuff::is_localizable_function<SourceImp>::value, "");
  static_assert(is_discrete_function<RangeImp>::value, "");
  typedef LocalizableOperatorDefault<GridViewImp, SourceImp, RangeImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::SourceType;
  using typename BaseType::RangeType;
  typedef typename Stuff::LA::Container<typename RangeType::RangeFieldType, container_backend>::MatrixType MatrixType;
  typedef typename Stuff::LA::Container<typename RangeType::RangeFieldType, container_backend>::VectorType VectorType;

private:
  typedef L2MatrixOperator<typename RangeType::SpaceType, MatrixType, GridViewType> LhsOperatorType;
  typedef L2VolumeVectorFunctional<SourceType, typename RangeType::SpaceType, VectorType, GridViewType>
      RhsFunctionalType;

public:
  L2GlobalProjectionLocalizableOperator(GridViewType grd_vw, const SourceType& src, RangeType& rng,
                                        const size_t over_integrate = 0)
    : BaseType(grd_vw, src, rng)
    , lhs_operator_(over_integrate, range_.space(), grid_view_)
    , rhs_functional_(over_integrate, source_, range_.space(), grid_view_)
    , solved_(false)
  {
    this->add(lhs_operator_);
    this->add(rhs_functional_);
    issue_warning(this->range().space());
  }

  void apply()
  {
    if (solved_)
      return;
    BaseType::apply();
    try {
      Stuff::LA::Solver<MatrixType>(lhs_operator_.matrix()).apply(rhs_functional_.vector(), range_.vector());
    } catch (Stuff::Exceptions::linear_solver_failed& ee) {
      DUNE_THROW(Exceptions::projection_error,
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
      DSC::TimedLogger().get("gdt.l2globalprojectionlocalizableoperator").warn()
          << "You are using this operator to project onto a discontinuous discrete function space (see below)!\n"
          << "Consider to use L2LocalProjectionLocalizableOperator instead!\n"
          << "You can disable this warning by defining "
          << "DUNE_GDT_OPERATORS_PROJECTIONS_L2_GLOBAL_LOCALIZABLE_DISABLE_WARNING\n"
          << "at compile time or by disabling the Dune::Stuff::Common::TimedLogger() instance at runtime.\n"
          << "The type of the range space is: " << DSC::Typename<S>::value() << std::endl;
    } // ... issue_warning(...)
  };

  template <class S>
  static inline void issue_warning(const S&)
  {
#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_L2_GLOBAL_LOCALIZABLE_DISABLE_WARNING
    Warning<S, !S::continuous>::issue();
#endif
  }

  using BaseType::grid_view_;
  using BaseType::source_;
  using BaseType::range_;

  LhsOperatorType lhs_operator_;
  RhsFunctionalType rhs_functional_;
  bool solved_;
}; // class L2GlobalProjectionLocalizableOperator


template <class GridViewType, class SourceType, class SpaceType, class VectorType>
typename std::enable_if<Stuff::Grid::is_grid_layer<GridViewType>::value
                            && Stuff::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && Stuff::LA::is_vector<VectorType>::value,
                        std::unique_ptr<L2GlobalProjectionLocalizableOperator<GridViewType, SourceType,
                                                                              DiscreteFunction<SpaceType, VectorType>,
                                                                              VectorType::sparse_matrix_type>>>::type
make_global_l2_projection_localizable_operator(const GridViewType& grid_view, const SourceType& source,
                                               DiscreteFunction<SpaceType, VectorType>& range,
                                               const size_t over_integrate = 0)
{
  return DSC::make_unique<L2GlobalProjectionLocalizableOperator<GridViewType,
                                                                SourceType,
                                                                DiscreteFunction<SpaceType, VectorType>,
                                                                VectorType::sparse_matrix_type>>(
      grid_view, source, range, over_integrate);
} // ... make_global_l2_projection_localizable_operator(...)

template <class SourceType, class SpaceType, class VectorType>
typename std::enable_if<Stuff::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && Stuff::LA::is_vector<VectorType>::value,
                        std::unique_ptr<L2GlobalProjectionLocalizableOperator<
                            typename SpaceType::GridViewType, SourceType, DiscreteFunction<SpaceType, VectorType>,
                            VectorType::sparse_matrix_type>>>::type
make_global_l2_projection_localizable_operator(const SourceType& source, DiscreteFunction<SpaceType, VectorType>& range,
                                               const size_t over_integrate = 0)
{
  return DSC::make_unique<L2GlobalProjectionLocalizableOperator<typename SpaceType::GridViewType,
                                                                SourceType,
                                                                DiscreteFunction<SpaceType, VectorType>,
                                                                VectorType::sparse_matrix_type>>(
      range.space().grid_view(), source, range, over_integrate);
} // ... make_global_l2_projection_localizable_operator(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROJECTIONS_L2_GLOBAL_HH
