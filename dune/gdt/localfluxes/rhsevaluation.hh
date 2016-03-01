#ifndef DUNE_GDT_LOCALFLUXES_RHSEVALUATION_HH
#define DUNE_GDT_LOCALFLUXES_RHSEVALUATION_HH

#include <dune/gdt/localevaluation/interface.hh>

namespace Dune {
namespace GDT {


// forward
template< class RHSFunctionImp >
class RHSEvaluation;


namespace internal {


template< class RHSFunctionImp >
class RHSEvaluationTraits
{
public:
  typedef RHSFunctionImp                                            RHSFunctionType;
  typedef RHSEvaluation< RHSFunctionType >                       derived_type;
  typedef typename RHSFunctionType::EntityType                      EntityType;
  typedef typename RHSFunctionType::LocalfunctionType               LocalfunctionType;
  typedef std::tuple< typename std::unique_ptr< LocalfunctionType > >  LocalfunctionTupleType;
  typedef typename RHSFunctionType::DomainFieldType                 DomainFieldType;
  static const size_t dimDomain = RHSFunctionType::dimDomain;
  static const size_t dimRange = RHSFunctionType::dimRange;
};


} // namespace internal


template< class RHSFunctionImp >
class RHSEvaluation
  : public Local
{
public:
  typedef typename internal::RHSEvaluationTraits< RHSFunctionImp > Traits;
  typedef typename Traits::RHSFunctionType                            RHSFunctionType;
  typedef typename Traits::EntityType                                    EntityType;
  typedef typename Traits::LocalfunctionTupleType                        LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType                               DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit RHSEvaluation(const RHSFunctionType& rhs_function)
    : rhs_function_(rhs_function)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(rhs_function_.local_global_function(entity));
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   */
  template< class R, size_t r, size_t rC >
  size_t order(const LocalfunctionTupleType& /*localFunctions_in*/,
               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& /*testBase*/)
  const
  {
    DUNE_THROW(NotImplemented, "Not meant to be integrated");
  }

  /**
   *  \brief  Computes a unary codim 0 evaluation.
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeCols of the testBase
   *  \attention ret is assumed to be zero!
   */
  template< class R, size_t rC >
  void evaluate(const LocalfunctionTupleType& local_source_function,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, rC >& entityAverage,
                const Dune::FieldVector< DomainFieldType, dimDomain >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
#if DUNE_VERSION_NEWER(DUNE_COMMON,3,9) //EXADUNE
    ret = Dune::DynamicVector< R >(std::get< 0 >(local_source_function)->evaluate(localPoint, entityAverage.evaluate(localPoint)[0]));
#else
    const auto fieldvector_ret = std::get< 0 >(local_source_function)->evaluate(localPoint,
                                                                                entityAverage.evaluate(localPoint)[0]);
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = fieldvector_ret[ii];
#endif
  }

private:
  const SourceFunctionType& source_function_;
}; // class Codim0Interface< Traits, 1 >


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFLUXES_RHSEVALUATION_HH
