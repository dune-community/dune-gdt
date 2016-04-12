#ifndef DUNE_GDT_LOCAL_FLUXES_INTERFACES_HH
#define DUNE_GDT_LOCAL_FLUXES_INTERFACES_HH

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace GDT {

namespace internal {


class IsNumericalCouplingFlux
{
};

class IsNumericalBoundaryFlux
{
};

class IsAnalyticalFlux
{
};

class IsRHSEvaluation
{
};


} // namespace internal


template <class Traits>
class LocalNumericalCouplingFluxInterface
    : public Stuff::CRTPInterface<LocalNumericalCouplingFluxInterface<Traits>, Traits>,
      internal::IsNumericalCouplingFlux
{
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;

public:
  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().local_functions(entity))
    return this->as_imp().local_functions(entity);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType>
  auto evaluate(const LocalfunctionTupleType& local_functions_tuple_entity,
                const LocalfunctionTupleType& local_functions_tuple_neighbor,
                const Stuff::LocalfunctionInterface<E, D, d, R, r, rC>& local_source_entity,
                const Stuff::LocalfunctionInterface<E, D, d, R, r, rC>& local_source_neighbor,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& x_intersection) const ->
      typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType
  {
    CHECK_CRTP(this->as_imp().evaluate(local_functions_tuple_entity,
                                       local_functions_tuple_neighbor,
                                       local_source_entity,
                                       local_source_neighbor,
                                       intersection,
                                       x_intersection));
    this->as_imp().evaluate(local_functions_tuple_entity,
                            local_functions_tuple_neighbor,
                            local_source_entity,
                            local_source_neighbor,
                            intersection,
                            x_intersection);
  }
}; // class LocalNumericalCouplingFluxInterface


template <class Traits>
class LocalNumericalBoundaryFluxInterface
    : public Stuff::CRTPInterface<LocalNumericalBoundaryFluxInterface<Traits>, Traits>,
      internal::IsNumericalBoundaryFlux
{
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;

public:
  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().local_functions(entity));
    return this->as_imp().local_functions(entity);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC, class IntersectionType>
  auto evaluate(const LocalfunctionTupleType& local_functions_tuple,
                const Stuff::LocalfunctionInterface<E, D, d, R, r, rC>& local_source_entity,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& x_intersection) const ->
      typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType
  {
    CHECK_CRTP(this->as_imp().evaluate(local_functions_tuple, local_source_entity, intersection, x_intersection))
    return this->as_imp().evaluate(local_functions_tuple, local_source_entity, intersection, x_intersection);
  }
}; // class LocalNumericalBoundaryFluxInterface


/** Analytical flux f for problem of the form delta_t u + div f(u,x,t) = 0 where u: R^d \to R^{r \times rC}.
 *  TODO: implement for rC > 1.
 *  TODO: determine correct FluxJacobianRangeType. */
template <class E, class D, size_t d, class R, size_t r, size_t rC = 1, bool gradient = false>
class AnalyticalFluxInterface : internal::IsAnalyticalFlux
{
  static_assert(rC == 1, "Not implemented for rC > 1");

public:
  typedef E EntityType;
  typedef D DomainFieldType;
  typedef R RangeFieldType;
  static const size_t dimDomain    = d;
  static const size_t dimRange     = r;
  static const size_t dimRangeCols = rC;

  virtual ~AnalyticalFluxInterface() = default;

  // Arbitrary entity type with dimension r for FluxRangeType and FluxJacobianRangeType definitions
  typedef typename Dune::template YaspGrid<r>::template Codim<0>::Entity FluxDummyEntityType;
  typedef typename Stuff::LocalfunctionSetInterface<FluxDummyEntityType, D, r, R, r, d>::RangeType FluxRangeType;
  // TODO: determine correct type
  typedef typename Stuff::LocalfunctionSetInterface<FluxDummyEntityType, D, r, R, r, d>::JacobianRangeType
      FluxJacobianRangeType;

  typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType RangeType; // of u, FieldVector or
  // FieldMatrix depending on
  // dimensions
  typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::DomainType DomainType;

  virtual FluxRangeType evaluate(const RangeType& u, const E& entity, const DomainType& x_local,
                                 const double t_ = 0) const = 0;

  virtual FluxJacobianRangeType jacobian(const RangeType& u, const E& entity, const DomainType& x_local,
                                         const double t_ = 0) const = 0;
}; // class AnalyticalFluxInterface<..., false>


/** Analytical flux f for problem of the form delta_t u + div f(u,x,t,\nabla u) = 0 where u: R^d \to R^r \times rC}.
 *  TODO: implement for rC > 1.
 *  TODO: determine correct FluxJacobianRangeType. */
template <class E, class D, size_t d, class R, size_t r, size_t rC>
class AnalyticalFluxInterface<E, D, d, R, r, rC, true> : internal::IsAnalyticalFlux
{
  static_assert(rC == 1, "Not implemented for rC > 1");

public:
  typedef E EntityType;
  typedef D DomainFieldType;
  typedef R RangeFieldType;
  static const size_t dimDomain    = d;
  static const size_t dimRange     = r;
  static const size_t dimRangeCols = rC;

  // Arbitrary entity type with dimension r for FluxRangeType and FluxJacobianRangeType definitions
  typedef typename Dune::template YaspGrid<r>::template Codim<0>::Entity FluxDummyEntityType;

  typedef typename Stuff::LocalfunctionSetInterface<FluxDummyEntityType, D, r, R, r, d>::RangeType FluxRangeType;
  // TODO: determine correct type
  typedef typename Stuff::LocalfunctionSetInterface<FluxDummyEntityType, D, r, R, r, d>::JacobianRangeType
      FluxJacobianRangeType;

  typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType RangeType; // of u, FieldVector or
  // FieldMatrix depending on
  // dimensions
  typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::JacobianRangeType
      JacobianRangeType; // of \nabla u, FieldMatrix or FieldVector< FieldMatrix > depending on dimensions
  typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::DomainType DomainType;

  virtual FluxRangeType evaluate(const RangeType& u, const JacobianRangeType& grad_u, const E& entity,
                                 const DomainType& x_local, const double t_ = 0) const = 0;

  virtual FluxJacobianRangeType jacobian(const RangeType& u, const JacobianRangeType& grad_u, const E& entity,
                                         const DomainType& x_local, const double t_ = 0) const = 0;
}; // class AnalyticalFluxInterface<..., true>


/**
 * Interface for analytical fluxes f that do not depend explicitly on the spatial variable x, the temporal variable t or
 * the gradient \nabla u. The corresponding PDE takes the form delta_t u + div f(u) = 0 where u: R^d \to
 * R^{r \times rC}. The flux f(u) is a function f: R^{r \times rC} \to R^{r \times rC \times d}.
 * TODO: implement for rC > 1.
 * */
template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class AutonomousAnalyticalFluxInterface : public AnalyticalFluxInterface<E, D, d, R, r, rC, false>
{
  static_assert(rC == 1, "Not implemented for rC > 1");
  typedef AnalyticalFluxInterface<E, D, d, R, r, rC, false> BaseType;

public:
  using typename BaseType::RangeType;
  using typename BaseType::FluxRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::FluxJacobianRangeType;

  virtual FluxRangeType evaluate(const RangeType& u, const E& /*entity*/, const DomainType& /*x_local*/,
                                 const double /*t_*/ = 0) const
  {
    return evaluate(u);
  }

  virtual FluxJacobianRangeType jacobian(const RangeType& u, const E& /*entity*/, const DomainType& /*x_local*/,
                                         const double /*t_*/ = 0) const
  {
    return jacobian(u);
  }

  virtual FluxRangeType evaluate(const RangeType& u) const = 0;
  virtual FluxJacobianRangeType jacobian(const RangeType& u) const = 0;
};


/**
 * Interface for right-hand side q of the PDE delta_t u + div f(u,x,t) = q(u,x,t) where u: R^d \to
 * R^{r \times rC}.
 * TODO: implement for rC > 1.
 * */
template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class RhsEvaluationFluxInterface : internal::IsRHSEvaluation
{
  static_assert(rC == 1, "Not implemented for rC > 1");

public:
  typedef E EntityType;
  typedef D DomainFieldType;
  typedef R RangeFieldType;
  static const size_t dimDomain    = d;
  static const size_t dimRange     = r;
  static const size_t dimRangeCols = rC;

  virtual ~RhsEvaluationFluxInterface() = default;

  typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType RangeType; // of u, FieldVector or
  // FieldMatrix depending on
  // dimensions
  typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>::DomainType DomainType;

  virtual RangeType evaluate(const RangeType& u, const E& entity, const DomainType& x_local,
                             const double t_ = 0) const = 0;
};


template <class T>
struct is_local_numerical_coupling_flux : std::is_base_of<internal::IsNumericalCouplingFlux, T>
{
};

template <class T>
struct is_local_numerical_boundary_flux : std::is_base_of<internal::IsNumericalBoundaryFlux, T>
{
};

template <class T>
struct is_analytical_flux : std::is_base_of<internal::IsAnalyticalFlux, T>
{
};

template <class T>
struct is_rhs_evaluation : std::is_base_of<internal::IsRHSEvaluation, T>
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_INTERFACES_HH
