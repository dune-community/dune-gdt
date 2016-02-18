#include <dune/gdt/localevaluation/interface.hh>


namespace Dune {
namespace GDT {


// forward
template< class SourceFunctionImp >
class SourceEvaluation;


namespace internal {


template< class SourceFunctionImp >
class SourceEvaluationTraits
{
public:
  typedef SourceFunctionImp                                            SourceFunctionType;
  typedef SourceEvaluation< SourceFunctionType >                       derived_type;
  typedef typename SourceFunctionType::EntityType                      EntityType;
  typedef typename SourceFunctionType::LocalfunctionType               LocalfunctionType;
  typedef std::tuple< typename std::unique_ptr< LocalfunctionType > >  LocalfunctionTupleType;
  typedef typename SourceFunctionType::DomainFieldType                 DomainFieldType;
  static const size_t dimDomain = SourceFunctionType::dimDomain;
  static const size_t dimRange = SourceFunctionType::rangeDimRange;
};


} // namespace internal


template< class SourceFunctionImp >
class SourceEvaluation
  : public LocalEvaluation::Codim0Interface< internal::SourceEvaluationTraits< SourceFunctionImp >, 1 >
{
public:
  typedef typename internal::SourceEvaluationTraits< SourceFunctionImp > Traits;
  typedef typename Traits::SourceFunctionType                            SourceFunctionType;
  typedef typename Traits::EntityType                                    EntityType;
  typedef typename Traits::LocalfunctionTupleType                        LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType                               DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit SourceEvaluation(const SourceFunctionType& source_function)
    : source_function_(source_function)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(source_function_.local_global_function(entity));
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
