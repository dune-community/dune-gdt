#ifndef DUNE_GDT_LOCALOPERATOR_FV_HH
#define DUNE_GDT_LOCALOPERATOR_FV_HH

#include "interfaces.hh"

namespace Dune {
namespace GDT {

// forwards
template <class NumericalFluxType>
class LocalCouplingFVOperator;

template <class NumericalFluxType>
class LocalBoundaryFVOperator;

template <class RHSEvaluationImp>
class LocalRHSFVOperator;


namespace internal {


// Traits
template <class NumericalFluxType>
struct LocalCouplingFVOperatorTraits
{
  typedef LocalCouplingFVOperator<NumericalFluxType> derived_type;
};

template <class NumericalFluxType>
struct LocalBoundaryFVOperatorTraits
{
  typedef LocalBoundaryFVOperator<NumericalFluxType> derived_type;
};

template <class RHSEvaluationImp>
struct LocalRHSFVOperatorTraits
{
  typedef LocalRHSFVOperator<RHSEvaluationImp> derived_type;
};


} // namespace internal


template <class NumericalFluxType>
class LocalCouplingFVOperator
    : public LocalCouplingOperatorInterface<internal::LocalCouplingFVOperatorTraits<NumericalFluxType>>
{
public:
  template <class... Args>
  explicit LocalCouplingFVOperator(Args&&... args)
    : numerical_flux_(std::forward<Args>(args)...)
  {
  }

  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source, const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor) const
  {
    const auto entity                         = intersection.inside();
    const auto neighbor                       = intersection.outside();
    const auto local_source_entity            = source.local_function(entity);
    const auto local_source_neighbor          = source.local_function(neighbor);
    const auto geometry_intersection          = intersection.geometry();
    const auto local_functions_tuple_entity   = numerical_flux_.local_functions(entity);
    const auto local_functions_tuple_neighbor = numerical_flux_.local_functions(neighbor);
    const auto result                         = numerical_flux_.evaluate(local_functions_tuple_entity,
                                                 local_functions_tuple_neighbor,
                                                 *local_source_entity,
                                                 *local_source_neighbor,
                                                 intersection,
                                                 geometry_intersection.local(geometry_intersection.center()));
    local_range_entity.vector().add(result / entity.geometry().volume());
    local_range_neighbor.vector().add(result * (-1.0 / neighbor.geometry().volume()));
  }

private:
  NumericalFluxType numerical_flux_;
};


template <class NumericalFluxType>
class LocalBoundaryFVOperator
    : public LocalBoundaryOperatorInterface<internal::LocalBoundaryFVOperatorTraits<NumericalFluxType>>
{
public:
  template <class... Args>
  explicit LocalBoundaryFVOperator(Args&&... args)
    : numerical_flux_(std::forward<Args>(args)...)
  {
  }

  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source, const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity) const
  {
    const auto entity                = intersection.inside();
    const auto local_source_entity   = source.local_function(entity);
    const auto geometry_intersection = intersection.geometry();
    const auto local_functions_tuple = numerical_flux_.local_functions(entity);
    const auto result                = numerical_flux_.evaluate(local_functions_tuple,
                                                 *local_source_entity,
                                                 intersection,
                                                 geometry_intersection.local(geometry_intersection.center()));
    local_range_entity.vector().add(result / entity.geometry().volume());
  }

private:
  NumericalFluxType numerical_flux_;
};

/** TODO: add support for time-dependent RHS
 *  TODO: implement as integral operator??
 * */
template <class RHSEvaluationImp>
class LocalRHSFVOperator : public LocalOperatorInterface<internal::LocalRHSFVOperatorTraits<RHSEvaluationImp>>
{
public:
  explicit LocalRHSFVOperator(const RHSEvaluationImp& rhs_evaluation)
    : rhs_evaluation_(rhs_evaluation)
  {
  }

  template <class SourceType, class RangeSpaceType, class VectorType>
  void apply(const SourceType& source, LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    const auto entity              = local_range.entity();
    const auto local_source_entity = source.local_function(entity);
    const auto x_local             = entity.geometry().local(entity.geometry().center());
    const auto u                   = local_source_entity->evaluate(x_local);
    const auto result = rhs_evaluation_.evaluate(u, entity, x_local);
    local_range.vector().add(result * entity.geometry().volume());
  }

private:
  const RHSEvaluationImp& rhs_evaluation_;
};

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_FV_HH
