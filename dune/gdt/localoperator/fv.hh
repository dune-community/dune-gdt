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


namespace internal {


// Traits
template <class NumericalFluxType>
class LocalCouplingFVOperatorTraits
{
  typedef LocalCouplingFVOperator<NumericalFluxType> derived_type;
};

template <class NumericalFluxType>
class LocalBoundaryFVOperatorTraits
{
  typedef LocalBoundaryFVOperator<NumericalFluxType> derived_type;
};


} // namespace internal


template <class NumericalFluxType>
class LocalBoundaryFVOperatorTraits;

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
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor)
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
                                                 local_source_entity,
                                                 local_source_neighbor,
                                                 intersection,
                                                 geometry_intersection.local(geometry_intersection.center()));
    local_range_entity += result;
    local_range_neighbor -= result;
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
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity)
  {
    const auto entity                = intersection.inside();
    const auto local_source_entity   = source.local_function(entity);
    const auto geometry_intersection = intersection.geometry();
    const auto local_functions_tuple = numerical_flux_.local_functions(entity);
    const auto result                = numerical_flux_.evaluate(local_functions_tuple,
                                                 local_source_entity,
                                                 intersection,
                                                 geometry_intersection.local(geometry_intersection.center()));
    local_range_entity += result;
  }

private:
  NumericalFluxType numerical_flux_;
};

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_FV_HH
