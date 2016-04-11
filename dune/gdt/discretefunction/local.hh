// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH
#define DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH

#include <vector>
#include <type_traits>

#include <dune/common/deprecated.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/fv/interface.hh>
#include <dune/gdt/mapper/interface.hh>

namespace Dune {
namespace GDT {


template <class SpaceImp, class VectorImp>
class ConstLocalDiscreteFunction
    : public Stuff::LocalfunctionInterface<typename SpaceImp::EntityType, typename SpaceImp::DomainFieldType,
                                           SpaceImp::dimDomain, typename SpaceImp::RangeFieldType, SpaceImp::dimRange,
                                           SpaceImp::dimRangeCols>
{
  static_assert(std::is_base_of<SpaceInterface<typename SpaceImp::Traits, SpaceImp::dimDomain, SpaceImp::dimRange,
                                               SpaceImp::dimRangeCols>,
                                SpaceImp>::value,
                "SpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<Dune::Stuff::LA::VectorInterface<typename VectorImp::Traits,
                                                                 typename VectorImp::Traits::ScalarType>,
                                VectorImp>::value,
                "VectorImp has to be derived from Stuff::LA::VectorInterface!");
  static_assert(std::is_same<typename SpaceImp::RangeFieldType, typename VectorImp::ScalarType>::value,
                "Types do not match!");
  typedef Stuff::LocalfunctionInterface<typename SpaceImp::EntityType, typename SpaceImp::DomainFieldType,
                                        SpaceImp::dimDomain, typename SpaceImp::RangeFieldType, SpaceImp::dimRange,
                                        SpaceImp::dimRangeCols> BaseType;
  typedef ConstLocalDiscreteFunction<SpaceImp, VectorImp> ThisType;

public:
  typedef SpaceImp SpaceType;
  typedef VectorImp VectorType;
  typedef typename BaseType::EntityType EntityType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const size_t dimRangeRows = BaseType::dimRangeCols;
  static const size_t dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;

private:
  typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;

public:
  typedef ConstLocalDoFVector<VectorType> ConstLocalDoFVectorType;

  ConstLocalDiscreteFunction(const SpaceType& sp, const VectorType& globalVector, const EntityType& ent)
    : BaseType(ent)
    , space_(sp)
    , base_(new BaseFunctionSetType(space_.base_function_set(this->entity())))
    , localVector_(new ConstLocalDoFVectorType(space_.mapper(), this->entity(), globalVector))
  {
    assert(localVector_->size() == base_->size());
  }

  ConstLocalDiscreteFunction(ThisType&& source) = default;

  ConstLocalDiscreteFunction(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  virtual ~ConstLocalDiscreteFunction()
  {
  }

  const SpaceType& space() const
  {
    return space_;
  }

  const BaseFunctionSetType& basis() const
  {
    return *base_;
  }

  const ConstLocalDoFVectorType& vector() const
  {
    return *localVector_;
  }

  virtual size_t order() const override
  {
    return base_->order();
  }

  void evaluate(const DomainType& xx, RangeType& ret) const override final
  {
    assert(this->is_a_valid_point(xx));
    if (!GDT::is_fv_space<SpaceType>::value) {
      std::fill(ret.begin(), ret.end(), RangeFieldType(0));
      std::vector<RangeType> tmpBaseValues(base_->size(), RangeType(0));
      assert(localVector_->size() == tmpBaseValues.size());
      base_->evaluate(xx, tmpBaseValues);
      for (size_t ii = 0; ii < localVector_->size(); ++ii) {
        ret.axpy(localVector_->get(ii), tmpBaseValues[ii]);
      }
    } else {
      for (size_t ii = 0; ii < localVector_->size(); ++ii)
        ret[ii] = localVector_->get(ii);
    }
  } // ... evaluate(...)

  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override final
  {
    assert(this->is_a_valid_point(xx));
    if (!GDT::is_fv_space<SpaceType>::value) {
      std::fill(ret.begin(), ret.end(), RangeFieldType(0));
      std::vector<JacobianRangeType> tmpBaseJacobianValues(base_->size(), JacobianRangeType(0));
      assert(localVector_->size() == tmpBaseJacobianValues.size());
      base_->jacobian(xx, tmpBaseJacobianValues);
      for (size_t ii = 0; ii < localVector_->size(); ++ii)
        ret.axpy(localVector_->get(ii), tmpBaseJacobianValues[ii]);
    } else {
      ret = JacobianRangeType(0);
    }
  } // ... jacobian(...)

  using BaseType::evaluate;
  using BaseType::jacobian;

protected:
  const SpaceType& space_;
  std::unique_ptr<const BaseFunctionSetType> base_;
  std::unique_ptr<const ConstLocalDoFVectorType> localVector_;
}; // class ConstLocalDiscreteFunction


template <class SpaceImp, class VectorImp>
class LocalDiscreteFunction : public ConstLocalDiscreteFunction<SpaceImp, VectorImp>
{
  typedef ConstLocalDiscreteFunction<SpaceImp, VectorImp> BaseType;
  typedef LocalDiscreteFunction<SpaceImp, VectorImp> ThisType;

public:
  typedef typename BaseType::SpaceType SpaceType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::EntityType EntityType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const size_t dimRangeRows = BaseType::dimRangeCols;
  static const size_t dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;

private:
  typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;

public:
  typedef LocalDoFVector<VectorType> LocalDoFVectorType;

  LocalDiscreteFunction(const SpaceType& sp, VectorType& globalVector, const EntityType& ent)
    : BaseType(sp, globalVector, ent)
    , localVector_(new LocalDoFVectorType(space_.mapper(), entity_, globalVector))
  {
    assert(localVector_->size() == base_->size());
  }

  //! previous comment questioned validity, defaulting this doesn't touch that question
  LocalDiscreteFunction(ThisType&& source) = default;

  LocalDiscreteFunction(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  virtual ~LocalDiscreteFunction()
  {
  }

  LocalDoFVectorType& vector()
  {
    return *localVector_;
  }

private:
  using BaseType::space_;
  using BaseType::entity_;
  using BaseType::base_;
  std::unique_ptr<LocalDoFVectorType> localVector_;
}; // class LocalDiscreteFunction


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH
