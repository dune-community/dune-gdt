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

#include "local-assemblers.hh"
#include "wrapper.hh"

namespace Dune {
namespace GDT {


template <class TestSpaceImp,
          class GridLayerImp = typename TestSpaceImp::GridLayerType,
          class AnsatzSpaceImp = TestSpaceImp,
          class OuterTestSpaceImp = TestSpaceImp,
          class OuterAnsatzSpaceImp = AnsatzSpaceImp>
class SystemAssembler : public XT::Grid::Walker<GridLayerImp>
{
  static_assert(is_space<TestSpaceImp>::value, "");
  static_assert(is_space<AnsatzSpaceImp>::value, "");
  static_assert(is_space<OuterTestSpaceImp>::value, "");
  static_assert(is_space<OuterAnsatzSpaceImp>::value, "");
  static_assert(std::is_same<typename TestSpaceImp::EntityType, XT::Grid::extract_entity_t<GridLayerImp>>::value, "");
  static_assert(std::is_same<typename AnsatzSpaceImp::EntityType, XT::Grid::extract_entity_t<GridLayerImp>>::value, "");
  static_assert(std::is_same<typename OuterTestSpaceImp::EntityType, XT::Grid::extract_entity_t<GridLayerImp>>::value,
                "");
  static_assert(std::is_same<typename OuterAnsatzSpaceImp::EntityType, XT::Grid::extract_entity_t<GridLayerImp>>::value,
                "");
  typedef XT::Grid::Walker<GridLayerImp> BaseType;
  typedef SystemAssembler<TestSpaceImp, GridLayerImp, AnsatzSpaceImp, OuterTestSpaceImp, OuterAnsatzSpaceImp> ThisType;

public:
  typedef TestSpaceImp TestSpaceType;
  typedef AnsatzSpaceImp AnsatzSpaceType;
  typedef OuterTestSpaceImp OuterTestSpaceType;
  typedef OuterAnsatzSpaceImp OuterAnsatzSpaceType;

  typedef typename BaseType::GridLayerType GridLayerType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  typedef XT::Grid::ApplyOn::WhichEntity<GridLayerType> ApplyOnWhichEntity;
  typedef XT::Grid::ApplyOn::WhichIntersection<GridLayerType> ApplyOnWhichIntersection;

  template <typename TestSpace,
            typename AnsatzSpace,
            typename = typename std::enable_if<std::is_same<OuterTestSpaceType, TestSpace>::value
                                               && std::is_same<OuterAnsatzSpaceType, AnsatzSpace>::value>::type>
  SystemAssembler(TestSpace test, AnsatzSpace ansatz, GridLayerType grd_layr)
    : BaseType(grd_layr)
    , test_space_(test)
    , ansatz_space_(ansatz)
    , outer_test_space_(test)
    , outer_ansatz_space_(ansatz)
  {
  }

  template <typename TestSpace,
            typename AnsatzSpace,
            typename =
                typename std::enable_if<std::is_same<OuterTestSpaceType, TestSpace>::value
                                        && std::is_same<OuterAnsatzSpaceType, AnsatzSpace>::value
                                        && std::is_same<typename TestSpace::GridLayerType, GridLayerType>::value>::type>
  explicit SystemAssembler(TestSpace test, AnsatzSpace ansatz)
    : BaseType(test.grid_layer())
    , test_space_(test)
    , ansatz_space_(ansatz)
    , outer_test_space_(test)
    , outer_ansatz_space_(ansatz)
  {
  }

  template <typename TestSpace,
            typename =
                typename std::enable_if<std::is_same<AnsatzSpaceType, TestSpace>::value
                                        && std::is_same<OuterTestSpaceType, TestSpace>::value
                                        && std::is_same<OuterAnsatzSpaceType, TestSpace>::value
                                        && std::is_same<typename TestSpace::GridLayerType, GridLayerType>::value>::type>
  explicit SystemAssembler(TestSpace test)
    : BaseType(test.grid_layer())
    , test_space_(test)
    , ansatz_space_(test)
    , outer_test_space_(test)
    , outer_ansatz_space_(test)
  {
  }

  template <typename TestSpace,
            typename = typename std::enable_if<std::is_same<AnsatzSpaceType, TestSpace>::value
                                               && std::is_same<OuterTestSpaceType, TestSpace>::value
                                               && std::is_same<OuterAnsatzSpaceType, TestSpace>::value>::type>
  explicit SystemAssembler(TestSpace test, GridLayerType grd_layr)
    : BaseType(grd_layr)
    , test_space_(test)
    , ansatz_space_(test)
    , outer_test_space_(test)
    , outer_ansatz_space_(test)
  {
  }

  SystemAssembler(GridLayerType grd_layr,
                  TestSpaceType inner_test,
                  AnsatzSpaceType inner_ansatz,
                  OuterTestSpaceType outer_test,
                  OuterAnsatzSpaceType outer_ansatz)
    : BaseType(grd_layr)
    , test_space_(inner_test)
    , ansatz_space_(inner_ansatz)
    , outer_test_space_(outer_test)
    , outer_ansatz_space_(outer_ansatz)
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

  const OuterTestSpaceType& outer_test_space() const
  {
    return *outer_test_space_;
  }

  const OuterAnsatzSpaceType& outer_ansatz_space() const
  {
    return *outer_ansatz_space_;
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

  template <class M, class R>
  ThisType& DUNE_DEPRECATED_MSG("Directly append the LocalVolumeTwoForm (13.05.2017)!") append(
      const LocalVolumeTwoFormAssembler<TestSpaceType, typename M::derived_type, AnsatzSpaceType>& local_assembler,
      XT::LA::MatrixInterface<M, R>& matrix,
      const ApplyOnWhichEntity* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalVolumeTwoFormMatrixAssemblerWrapper<ThisType, typename M::derived_type> WrapperType;
    this->codim0_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
    return *this;
  } // ... append(...)

  template <class M, class R>
  ThisType& append(const LocalVolumeTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                                     typename AnsatzSpaceType::BaseFunctionSetType,
                                                     typename M::ScalarType>& local_volume_two_form,
                   XT::LA::MatrixInterface<M, R>& matrix,
                   const ApplyOnWhichEntity* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    this->codim0_functors_.emplace_back(
        new LocalVolumeTwoFormAssemblerFunctor<TestSpaceType, typename M::derived_type, GridLayerType, AnsatzSpaceType>(
            test_space_, ansatz_space_, where, local_volume_two_form, matrix.as_imp()));
    return *this;
  } // ... append(...)

  template <class M, class R>
  ThisType& append(const LocalCouplingTwoFormAssembler<TestSpaceType,
                                                       IntersectionType,
                                                       typename M::derived_type,
                                                       AnsatzSpaceType,
                                                       OuterTestSpaceType,
                                                       OuterAnsatzSpaceType>& local_assembler,
                   XT::LA::MatrixInterface<M, R>& matrix,
                   const ApplyOnWhichIntersection* where = new XT::Grid::ApplyOn::AllIntersections<GridLayerType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalCouplingTwoFormMatrixAssemblerWrapper<ThisType, typename M::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
    return *this;
  } // ... append(...)

  template <class M, class R>
  ThisType& append(const LocalCouplingTwoFormAssembler<TestSpaceType,
                                                       IntersectionType,
                                                       typename M::derived_type,
                                                       AnsatzSpaceType,
                                                       OuterTestSpaceType,
                                                       OuterAnsatzSpaceType>& local_assembler,
                   XT::LA::MatrixInterface<M, R>& matrix_in_in,
                   XT::LA::MatrixInterface<M, R>& matrix_out_out,
                   XT::LA::MatrixInterface<M, R>& matrix_in_out,
                   XT::LA::MatrixInterface<M, R>& matrix_out_in,
                   const ApplyOnWhichIntersection* where = new XT::Grid::ApplyOn::AllIntersections<GridLayerType>())
  {
    assert(matrix_in_in.rows() == test_space_->mapper().size());
    assert(matrix_in_in.cols() == ansatz_space_->mapper().size());
    assert(matrix_out_out.rows() == outer_test_space_->mapper().size());
    assert(matrix_out_out.cols() == outer_ansatz_space_->mapper().size());
    assert(matrix_in_out.rows() == test_space_->mapper().size());
    assert(matrix_in_out.cols() == outer_ansatz_space_->mapper().size());
    assert(matrix_out_in.rows() == outer_test_space_->mapper().size());
    assert(matrix_out_in.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalCouplingTwoFormMatrixAssemblerWrapper<ThisType, typename M::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_,
                                                        ansatz_space_,
                                                        outer_test_space_,
                                                        outer_ansatz_space_,
                                                        where,
                                                        local_assembler,
                                                        matrix_in_in.as_imp(),
                                                        matrix_out_out.as_imp(),
                                                        matrix_in_out.as_imp(),
                                                        matrix_out_in.as_imp()));
    return *this;
  } // ... append(...)

  template <class M, class R>
  ThisType& append(
      const LocalBoundaryTwoFormAssembler<TestSpaceType, IntersectionType, typename M::derived_type, AnsatzSpaceType>&
          local_assembler,
      XT::LA::MatrixInterface<M, R>& matrix,
      const ApplyOnWhichIntersection* where = new XT::Grid::ApplyOn::AllIntersections<GridLayerType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalBoundaryTwoFormMatrixAssemblerWrapper<ThisType, typename M::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
    return *this;
  } // ... append(...)

  template <class V, class R>
  ThisType& append(const LocalVolumeFunctionalAssembler<TestSpaceType, typename V::derived_type>& local_assembler,
                   XT::LA::VectorInterface<V, R>& vector,
                   const ApplyOnWhichEntity* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalVolumeFunctionalVectorAssemblerWrapper<ThisType, typename V::derived_type> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
    return *this;
  } // ... append(...)

  template <class V, class R>
  ThisType&
  append(const LocalFaceFunctionalAssembler<TestSpaceType, IntersectionType, typename V::derived_type>& local_assembler,
         XT::LA::VectorInterface<V, R>& vector,
         const ApplyOnWhichIntersection* where = new XT::Grid::ApplyOn::AllIntersections<GridLayerType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalFaceFunctionalVectorAssemblerWrapper<ThisType, typename V::derived_type> WrapperType;
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
  const Dune::XT::Common::PerThreadValue<const OuterTestSpaceType> outer_test_space_;
  const Dune::XT::Common::PerThreadValue<const OuterAnsatzSpaceType> outer_ansatz_space_;
}; // class SystemAssembler


} // namespace GDT
} // namespace Dune


#include "system.lib.hh"


#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
