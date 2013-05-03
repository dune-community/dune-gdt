#ifndef DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONSTRAINTS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONSTRAINTS_HH

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/static_assert.hh>

#include "interface.hh"
#include "lagrange-continuous.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


template <class Traits>
class ConstraintsInterface;


template <class RangeFieldImp = double>
class LocalConstraints
{
public:
  typedef RangeFieldImp RangeFieldType;

  typedef Dune::DynamicVector<size_t> IndicesType;
  typedef Dune::DynamicMatrix<RangeFieldType> ValuesType;

public:
  LocalConstraints(const size_t numRows, const size_t numCols)
    : globalRows_(numRows)
    , globalCols_(numCols)
    , values_(numRows, numCols)
  {
  }

  size_t rows() const
  {
    return globalRows_.size();
  }

  size_t cols() const
  {
    return globalCols_.size();
  }

  const IndicesType& globalRows() const
  {
    return globalRows_;
  }

  const IndicesType& globalCols() const
  {
    return globalCols_;
  }

  const ValuesType& values() const
  {
    return values_;
  }

private:
  template <class Traits>
  friend class ConstraintsInterface;

  IndicesType globalRows_;
  IndicesType globalCols_;
  ValuesType values_;
}; // class LocalConstraints


template <class Traits>
class ConstraintsInterface
{
public:
  typedef typename Traits::derived_type derived_type;

  template <class T, class A>
  LocalConstraints<typename SpaceInterface<T>::RangeFieldType>
  local(const SpaceInterface<T>& testSpace, const SpaceInterface<A>& ansatzSpace,
        const typename SpaceInterface<T>::EntityType& entity) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().local(testSpace, ansatzSpace, entity));
    return asImp().local(testSpace, ansatzSpace, entity);
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class ConstraintsInterface


// forward, to be used in the traits
class NoConstraints;

class NoConstraintsTraits
{
public:
  typedef NoConstraints derived_type;
};

class NoConstraints : public ConstraintsInterface<NoConstraintsTraits>
{
public:
  typedef NoConstraintsTraits Traits;

  template <class T, class A>
  LocalConstraints<typename SpaceInterface<T>::RangeFieldType>
  local(const SpaceInterface<T>& /*testSpace*/, const SpaceInterface<A>& /*ansatzSpace*/,
        const typename SpaceInterface<T>::EntityType& /*entity*/) const
  {
    return LocalConstraints<typename SpaceInterface<T>::RangeFieldType>(0, 0);
  }
};


// forward, to be used in the traits
class DirichletConstraints;

class DirichletConstraintsTraits
{
public:
  typedef DirichletConstraints derived_type;
};

class DirichletConstraints
{
public:
  typedef DirichletConstraintsTraits Traits;

  template <class T, class A>
  LocalConstraints<typename SpaceInterface<T>::RangeFieldType>
  local(const SpaceInterface<T>& testSpace, const SpaceInterface<A>& ansatzSpace,
        const typename SpaceInterface<T>::EntityType& entity) const
  {
    dune_static_assert((Dune::AlwaysFalse<T>::value), "ERROR: not implemeneted for arbitrary test spaces!");
  }
}; // class DirichletConstraints


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONSTRAINTS_HH
