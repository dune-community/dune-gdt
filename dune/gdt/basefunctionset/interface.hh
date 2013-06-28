#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_INTERFACE_HH

#include <vector>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {


/**
 *  \brief Interface for matrix valued basis functions.
 *
 *  \note   see specialization for rangeDimCols = 1 for vector and scalar valued and basis functions.
 */
template <class Traits, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDim = 1>
class BaseFunctionSetInterface;


/**
 *  \brief Interface for scalar valued basis functions.
 */
template <class Traits, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class BaseFunctionSetInterface<Traits, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;

  const EntityType& entity() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().entity());
    return asImp().entity();
  }

  const BackendType& backend() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().backend());
    return asImp().backend();
  }

  size_t size() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
    return asImp().size();
  }

  size_t order() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order());
    return asImp().order();
  }

  void evaluate(const DomainType& x, std::vector<RangeType>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(x, ret));
  }

  void jacobian(const DomainType& x, std::vector<JacobianRangeType>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().jacobian(x, ret));
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class BaseFunctionSetInterface


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_INTERFACE_HH
