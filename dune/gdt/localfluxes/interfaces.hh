#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace GDT {


template< class Traits >
class NumericalCouplingFluxInterface
    : Stuff::CRTPInterface< NumericalCouplingFluxInterface< Traits >, Traits >
{
  typedef typename Traits::ResultType ResultType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType             EntityType;

public:
  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().local_functions(entity));
    return this->as_imp().local_functions(entity);
  }

  template< class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType >
  ResultType evaluate(const LocalfunctionTupleType& local_functions_tuple_entity,
                      const LocalfunctionTupleType& local_functions_tuple_neighbor,
                      const Stuff::LocalfunctionInterface< E, D, d, R, r, rC >& local_source_entity,
                      const Stuff::LocalfunctionInterface< E, D, d, R, r, rC >& local_source_neighbor,
                      const IntersectionType& intersection,
                      const Dune::FieldVector< D, d - 1 >& x_intersection) const
  {
    CHECK_CRTP(this->as_imp().evaluate(local_functions_tuple_entity, local_functions_tuple_neigbhor, local_source_entity, local_source_neighbor, intersection, x_intersection));
    this->as_imp().evaluate(local_functions_tuple_entity, local_functions_tuple_neighbor, local_source_entity, local_source_neighbor, intersection, x_intersection);
  }
}; // class NumericalCouplingFluxInterface

template< class Traits >
class NumericalBoundaryFluxInterface
    : Stuff::CRTPInterface< NumericalBoundaryFluxInterface< Traits >, Traits >
{
  typedef typename Traits::ResultType ResultType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType             EntityType;
public:
  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().local_functions(entity));
    return this->as_imp().local_functions(entity);
  }

  template< class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType >
  ResultType evaluate(const LocalfunctionTupleType& local_functions_tuple,
                      const Stuff::LocalfunctionInterface< E, D, d, R, r, rC >& local_source_entity,
                      const IntersectionType& intersection,
                      const Dune::FieldVector< D, d - 1 >& x_intersection) const
  {
    CHECK_CRTP(this->as_imp().evaluate(local_functions_tuple, local_source_entity, intersection, x_intersection));
    this->as_imp().evaluate(local_functions_tuple, local_source_entity, local_source_neighbor, intersection, x_intersection);
  }
}; // class NumericalBoundaryFluxInterface


template< bool gradient = false >
class IsAnalyticalFlux
{
};

template< class E, class D, size_t d, class R, size_t r, size_t rC = 1, bool gradient = false>
class AnalyticalFluxInterface
    : IsAnalyticalFlux<gradient>
{
public:
  typedef ?? FluxRangeType;
  typedef ?? FluxJacobianRangeType;

  typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC>::RangeType RangeType; // of u, FieldVector or FieldMatrix depending on dimensions
  typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC>::DomainType DomainType;

  virtual FluxRangeType evaluate(const RangeType& u,
                  const E& entity,
                  const DomainType& x_local,
                  const double t_ = 0) const = 0;
}; // class AnalyticalFluxInterface<..., false>

template< class E, class D, size_t d, class R, size_t r, size_t rC >
class AnalyticalFluxInterface< E, D, d, R, r, rC, true >
    : IsAnalyticalFlux<true>
{
public:
  typedef E EntityType;
  typedef D DomainFieldType;
  typedef R RangeFieldType;

  typedef ?? FluxRangeType;
  typedef ?? FluxJacobianRangeType;

  typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC>::RangeType RangeType; // of u, FieldVector or FieldMatrix depending on dimensions
  typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC>::JacobianRangeType JacobianRangeType; // of gradient of u, FieldMatrix or FieldVector< FieldMatrix > depending on dimensions
  typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC>::DomainType DomainType;
  virtual FluxRangeType evaluate(const RangeType& u,
                  const JacobianRangeType& grad_u,
                  const E& entity,
                  const DomainType& x_local,
                  const double t_ = 0) const = 0;
}; // class AnalyticalFluxInterface<..., true>

class IsAutonomousAnalyticalFlux
    : IsAnalyticalFlux
{

};

/**
 * Interface for analytical fluxes that do not depend explicitly on the spatial variable x or the temporal variable t
 * */
template< class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class AutonomousAnalyticalFluxInterface
    : AnalyticalFluxInterface< E, D, d, R, r, rC, false >
    , IsAutonomousAnalyticalFlux
{
   typedef AnalyticalFluxInterface< E, D, d, R, r, rC, false > BaseType;
public:
  using BaseType::RangeType;
  using BaseType::FluxRangeType;
  using BaseType::DomainType;

  virtual FluxRangeType evaluate(const RangeType& u,
                  const E& /*entity*/,
                  const DomainType& /*x_local*/,
                  const double /*t_*/ = 0) const
  {
    return evaluate(u);
  }

  virtual FluxRangeType evaluate(const RangeType& u) const = 0;
};


} // namespace GDT
} // namespace Dune
