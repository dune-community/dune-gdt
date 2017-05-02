// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   Rene Milk       (2013 - 2014, 2016 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_HH

#include <type_traits>
#include <memory>

#include <dune/common/deprecated.hh>
#include <dune/common/version.hh>

#include <dune/xt/common/parallel/helper.hh>

#include <dune/xt/grid/walker.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>

#include "wrapper.hh"

namespace Dune {
namespace GDT {


template <class TestSpaceImp,
          class GridLayerImp = typename TestSpaceImp::GridLayerType,
          class AnsatzSpaceImp = TestSpaceImp>
class SystemAssembler : public XT::Grid::Walker<GridLayerImp>
{
  static_assert(GDT::is_space<TestSpaceImp>::value, "TestSpaceImp has to be derived from SpaceInterface!");
  static_assert(GDT::is_space<AnsatzSpaceImp>::value, "AnsatzSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_same<typename TestSpaceImp::RangeFieldType, typename AnsatzSpaceImp::RangeFieldType>::value,
                "Types do not match!");
  typedef XT::Grid::Walker<GridLayerImp> BaseType;
  typedef SystemAssembler<TestSpaceImp, GridLayerImp, AnsatzSpaceImp> ThisType;

public:
  typedef TestSpaceImp TestSpaceType;
  typedef AnsatzSpaceImp AnsatzSpaceType;
  typedef typename TestSpaceType::RangeFieldType RangeFieldType;

  typedef typename BaseType::GridLayerType GridLayerType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  typedef XT::Grid::ApplyOn::WhichEntity<GridLayerType> ApplyOnWhichEntity;
  typedef XT::Grid::ApplyOn::WhichIntersection<GridLayerType> ApplyOnWhichIntersection;

  SystemAssembler(TestSpaceType test, AnsatzSpaceType ansatz, GridLayerType grd_layr)
    : BaseType(grd_layr)
    , test_space_(test)
    , ansatz_space_(ansatz)
  {
  }

  /// \todo Guard against GridLayerType != TestSpaceImp::GridLayerType
  SystemAssembler(TestSpaceType test, AnsatzSpaceType ansatz)
    : BaseType(test.grid_layer())
    , test_space_(test)
    , ansatz_space_(ansatz)
  {
  }

  /// \todo Guard against AnsatzSpaceType != GridLayerType || GridLayerType != TestSpaceType::GridLayerType
  explicit SystemAssembler(TestSpaceType test)
    : BaseType(test.grid_layer())
    , test_space_(test)
    , ansatz_space_(test)
  {
  }

  /// \todo Guard against AnsatzSpaceType != TestSpaceType
  SystemAssembler(TestSpaceType test, GridLayerType grd_layr)
    : BaseType(grd_layr)
    , test_space_(test)
    , ansatz_space_(test)
  {
  }

  /// \sa https://github.com/dune-community/dune-gdt/issues/89
  SystemAssembler(const ThisType& other) = delete; // all wrappers hold references to dead spaces after move!
  SystemAssembler(ThisType&& source) = delete;
  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const TestSpaceType& test_space() const
  {
    return *test_space_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return *ansatz_space_;
  }

  using BaseType::append;

  template <class C>
  ThisType& append(ConstraintsInterface<C>& constraints,
                   const ApplyOnWhichEntity* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    typedef internal::ConstraintsWrapper<TestSpaceType, AnsatzSpaceType, GridLayerType, typename C::derived_type>
        WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, constraints.as_imp()));
    return *this;
  } // ... append(...)

  template <class V, class M>
  ThisType& append(const LocalVolumeTwoFormAssembler<V>& local_assembler,
                   XT::LA::MatrixInterface<M, RangeFieldType>& matrix,
                   const ApplyOnWhichEntity* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalVolumeTwoFormMatrixAssemblerWrapper<ThisType,
                                                               LocalVolumeTwoFormAssembler<V>,
                                                               typename M::derived_type>
        WrapperType;
    this->codim0_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
    return *this;
  } // ... append(...)

  template <class V, class R>
  ThisType& append(const LocalVolumeTwoFormAssembler<V>& local_assembler,
                   std::vector<DynamicMatrix<R>>& matrix,
                   const ApplyOnWhichEntity* where = new XT::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    assert(matrix.size() == test_space_->grid_view().size(0));
    typedef internal::LocalVolumeTwoFormMatrixAssemblerWrapper<ThisType,
                                                               LocalVolumeTwoFormAssembler<V>,
                                                               std::vector<DynamicMatrix<R>>>
        WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix));
    return *this;
  } // ... append(...)

  template <class V, class M>
  ThisType& append(const LocalCouplingTwoFormAssembler<V>& local_assembler,
                   XT::LA::MatrixInterface<M, RangeFieldType>& matrix,
                   const ApplyOnWhichIntersection* where = new XT::Grid::ApplyOn::AllIntersections<GridLayerType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalCouplingTwoFormMatrixAssemblerWrapper<ThisType,
                                                                 LocalCouplingTwoFormAssembler<V>,
                                                                 typename M::derived_type>
        WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
    return *this;
  } // ... append(...)

  template <class V, class M>
  ThisType& append(const LocalBoundaryTwoFormAssembler<V>& local_assembler,
                   XT::LA::MatrixInterface<M, RangeFieldType>& matrix,
                   const ApplyOnWhichIntersection* where = new XT::Grid::ApplyOn::AllIntersections<GridLayerType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalBoundaryTwoFormMatrixAssemblerWrapper<ThisType,
                                                                 LocalBoundaryTwoFormAssembler<V>,
                                                                 typename M::derived_type>
        WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
    return *this;
  } // ... append(...)

  template <class L, class V>
  ThisType& append(const LocalVolumeFunctionalAssembler<L>& local_assembler,
                   XT::LA::VectorInterface<V, RangeFieldType>& vector,
                   const ApplyOnWhichEntity* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalVolumeFunctionalVectorAssemblerWrapper<ThisType,
                                                                  LocalVolumeFunctionalAssembler<L>,
                                                                  typename V::derived_type>
        WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
    return *this;
  } // ... append(...)

  template <class L, class V>
  ThisType& append(const LocalFaceFunctionalAssembler<L>& local_assembler,
                   XT::LA::VectorInterface<V, RangeFieldType>& vector,
                   const ApplyOnWhichIntersection* where = new XT::Grid::ApplyOn::AllIntersections<GridLayerType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalFaceFunctionalVectorAssemblerWrapper<ThisType,
                                                                LocalFaceFunctionalAssembler<L>,
                                                                typename V::derived_type>
        WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
    return *this;
  } // ... append(...)

  using BaseType::add;

  template <class... Args>
  DUNE_DEPRECATED_MSG("Use append() instead (since 11.01.2017)!")
  ThisType& add(Args&&... args)
  {
    return append(std::forward<Args>(args)...);
  }

  void assemble(const bool use_tbb = false)
  {
    this->walk(use_tbb);
  }

  template <class Partitioning>
  void assemble(const Partitioning& partitioning)
  {
    this->walk(partitioning);
  }

protected:
  const Dune::XT::Common::PerThreadValue<const TestSpaceType> test_space_;
  const Dune::XT::Common::PerThreadValue<const AnsatzSpaceType> ansatz_space_;
}; // class SystemAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
