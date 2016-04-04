// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_HH

#include <type_traits>
#include <memory>

#include <dune/common/deprecated.hh>
#include <dune/common/version.hh>

#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/common/parallel/helper.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>

#include "local/codim0.hh"
#include "local/codim1.hh"
#include "wrapper.hh"

namespace Dune {
namespace GDT {


template <class TestSpaceImp, class GridViewImp = typename TestSpaceImp::GridViewType,
          class AnsatzSpaceImp                  = TestSpaceImp>
class SystemAssembler : public DSG::Walker<GridViewImp>
{
  static_assert(GDT::is_space<TestSpaceImp>::value, "TestSpaceImp has to be derived from SpaceInterface!");
  static_assert(GDT::is_space<AnsatzSpaceImp>::value, "AnsatzSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_same<typename TestSpaceImp::RangeFieldType, typename AnsatzSpaceImp::RangeFieldType>::value,
                "Types do not match!");
  typedef DSG::Walker<GridViewImp> BaseType;
  typedef SystemAssembler<TestSpaceImp, GridViewImp, AnsatzSpaceImp> ThisType;

public:
  typedef TestSpaceImp TestSpaceType;
  typedef AnsatzSpaceImp AnsatzSpaceType;
  typedef typename TestSpaceType::RangeFieldType RangeFieldType;

  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  typedef DSG::ApplyOn::WhichEntity<GridViewType> ApplyOnWhichEntity;
  typedef DSG::ApplyOn::WhichIntersection<GridViewType> ApplyOnWhichIntersection;

  SystemAssembler(TestSpaceType test, AnsatzSpaceType ansatz, GridViewType grd_vw)
    : BaseType(grd_vw)
    , test_space_(test)
    , ansatz_space_(ansatz)
  {
  }

  SystemAssembler(TestSpaceType test, AnsatzSpaceType ansatz)
    : BaseType(test.grid_view())
    , test_space_(test)
    , ansatz_space_(ansatz)
  {
  }

  explicit SystemAssembler(TestSpaceType test)
    : BaseType(test.grid_view())
    , test_space_(test)
    , ansatz_space_(test)
  {
  }

  SystemAssembler(TestSpaceType test, GridViewType grd_vw)
    : BaseType(grd_vw)
    , test_space_(test)
    , ansatz_space_(test)
  {
  }

  SystemAssembler(ThisType&& source) = default;

  const TestSpaceType& test_space() const
  {
    return *test_space_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return *ansatz_space_;
  }

  using BaseType::add;

  template <class C>
  void add(Spaces::ConstraintsInterface<C>& constraints,
           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    typedef internal::ConstraintsWrapper<TestSpaceType, AnsatzSpaceType, GridViewType, typename C::derived_type>
        WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, constraints.as_imp()));
  } // ... add(...)

  template <class L, class M>
  void add(const LocalAssembler::Codim0Matrix<L>& local_assembler,
           Stuff::LA::MatrixInterface<M, RangeFieldType>& matrix,
           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalVolumeMatrixAssemblerWrapper<ThisType,
                                                        LocalAssembler::Codim0Matrix<L>,
                                                        typename M::derived_type> WrapperType;
    this->codim0_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
  } // ... add(...)

  template <class V, class M>
  void add(const LocalVolumeTwoFormAssembler<V>& local_assembler, Stuff::LA::MatrixInterface<M, RangeFieldType>& matrix,
           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalVolumeTwoFormMatrixAssemblerWrapper<ThisType,
                                                               LocalVolumeTwoFormAssembler<V>,
                                                               typename M::derived_type> WrapperType;
    this->codim0_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
  } // ... add(...)

  template <class V, class M>
  void add(const LocalCouplingTwoFormAssembler<V>& local_assembler,
           Stuff::LA::MatrixInterface<M, RangeFieldType>& matrix,
           const ApplyOnWhichIntersection* where = new DSG::ApplyOn::AllIntersections<GridViewType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalCouplingTwoFormMatrixAssemblerWrapper<ThisType,
                                                                 LocalCouplingTwoFormAssembler<V>,
                                                                 typename M::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
  } // ... add(...)

  template <class V, class M>
  void add(const LocalBoundaryTwoFormAssembler<V>& local_assembler,
           Stuff::LA::MatrixInterface<M, RangeFieldType>& matrix,
           const ApplyOnWhichIntersection* where = new DSG::ApplyOn::AllIntersections<GridViewType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalBoundaryTwoFormMatrixAssemblerWrapper<ThisType,
                                                                 LocalBoundaryTwoFormAssembler<V>,
                                                                 typename M::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
  } // ... add(...)

  template <class L, class V>
  void add(const LocalVolumeFunctionalAssembler<L>& local_assembler,
           Stuff::LA::VectorInterface<V, RangeFieldType>& vector,
           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalVolumeFunctionalVectorAssemblerWrapper<ThisType,
                                                                  LocalVolumeFunctionalAssembler<L>,
                                                                  typename V::derived_type> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
  } // ... add(...)

  template <class L, class V>
  void add(const LocalFaceFunctionalAssembler<L>& local_assembler,
           Stuff::LA::VectorInterface<V, RangeFieldType>& vector,
           const ApplyOnWhichIntersection* where = new DSG::ApplyOn::AllIntersections<GridViewType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalFaceFunctionalVectorAssemblerWrapper<ThisType,
                                                                LocalFaceFunctionalAssembler<L>,
                                                                typename V::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
  } // ... add(...)

  template <class Codim0Assembler, class M>
  void DUNE_DEPRECATED_MSG("Will be removed or first argument has to be replaced by an interface (04.02.2015)!")
      add_codim0_assembler(const Codim0Assembler& local_assembler,
                           Stuff::LA::MatrixInterface<M, RangeFieldType>& matrix,
                           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalVolumeMatrixAssemblerWrapper<ThisType, Codim0Assembler, typename M::derived_type>
        WrapperType;
    this->codim0_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
  } // ... add(...)

  template <class Codim0Assembler, class V>
  void DUNE_DEPRECATED_MSG("Will be removed or first argument has to be replaced by an interface (04.02.2015)!")
      add_codim0_assembler(const Codim0Assembler& local_assembler,
                           Stuff::LA::VectorInterface<V, RangeFieldType>& vector,
                           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalVolumeVectorAssemblerWrapper<ThisType, Codim0Assembler, typename V::derived_type>
        WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
  } // ... add(...)

  template <class L, class M>
  void add(const LocalAssembler::Codim1CouplingMatrix<L>& local_assembler,
           Stuff::LA::MatrixInterface<M, RangeFieldType>& matrix,
           const ApplyOnWhichIntersection* where = new DSG::ApplyOn::AllIntersections<GridViewType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalFaceMatrixAssemblerWrapper<ThisType,
                                                      LocalAssembler::Codim1CouplingMatrix<L>,
                                                      typename M::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
  } // ... add(...)

  template <class L, class M>
  void add(const LocalAssembler::Codim1BoundaryMatrix<L>& local_assembler,
           Stuff::LA::MatrixInterface<M, RangeFieldType>& matrix,
           const ApplyOnWhichIntersection* where = new DSG::ApplyOn::AllIntersections<GridViewType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalFaceMatrixAssemblerWrapper<ThisType,
                                                      LocalAssembler::Codim1BoundaryMatrix<L>,
                                                      typename M::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
  } // ... add(...)

  template <class L, class V>
  void add(const LocalAssembler::Codim0Vector<L>& local_assembler,
           Stuff::LA::VectorInterface<V, RangeFieldType>& vector,
           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalVolumeVectorAssemblerWrapper<ThisType,
                                                        LocalAssembler::Codim0Vector<L>,
                                                        typename V::derived_type> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
  } // ... add(...)

  template <class L, class V>
  void add(const LocalAssembler::Codim1Vector<L>& local_assembler,
           Stuff::LA::VectorInterface<V, RangeFieldType>& vector,
           const ApplyOnWhichIntersection* where = new DSG::ApplyOn::AllIntersections<GridViewType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalFaceVectorAssemblerWrapper<ThisType,
                                                      LocalAssembler::Codim1Vector<L>,
                                                      typename V::derived_type> WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
  } // ... add(...)

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
  const DS::PerThreadValue<const TestSpaceType> test_space_;
  const DS::PerThreadValue<const AnsatzSpaceType> ansatz_space_;
}; // class SystemAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
