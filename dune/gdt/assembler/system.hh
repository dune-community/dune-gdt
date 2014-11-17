// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_HH

#include <type_traits>
#include <memory>

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


template< class TestSpaceImp,
          class GridViewImp = typename TestSpaceImp::GridViewType,
          class AnsatzSpaceImp = TestSpaceImp >
class SystemAssembler
  : public DSG::Walker< GridViewImp >
{
  static_assert(std::is_base_of< SpaceInterface< typename TestSpaceImp::Traits >, TestSpaceImp >::value,
                "TestSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of< SpaceInterface< typename AnsatzSpaceImp::Traits >, AnsatzSpaceImp >::value,
                "AnsatzSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_same< typename TestSpaceImp::RangeFieldType, typename AnsatzSpaceImp::RangeFieldType >::value,
                "Types do not match!");
  typedef DSG::Walker< GridViewImp >                                  BaseType;
  typedef SystemAssembler<TestSpaceImp, GridViewImp, AnsatzSpaceImp > ThisType;
public:
  typedef TestSpaceImp                           TestSpaceType;
  typedef AnsatzSpaceImp                         AnsatzSpaceType;
  typedef typename TestSpaceType::RangeFieldType RangeFieldType;

  typedef typename BaseType::GridViewType     GridViewType;
  typedef typename BaseType::EntityType       EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  typedef DSG::ApplyOn::WhichEntity< GridViewType >       ApplyOnWhichEntity;
  typedef DSG::ApplyOn::WhichIntersection< GridViewType > ApplyOnWhichIntersection;

  SystemAssembler(TestSpaceType test, AnsatzSpaceType ansatz, GridViewType grid_view)
    : BaseType(grid_view)
    , test_space_(test)
    , ansatz_space_(ansatz)
  {}

  SystemAssembler(TestSpaceType test, AnsatzSpaceType ansatz)
    : BaseType(test.grid_view())
    , test_space_(test)
    , ansatz_space_(ansatz)
  {}

  explicit SystemAssembler(TestSpaceType test)
    : BaseType(test.grid_view())
    , test_space_(test)
    , ansatz_space_(test)
  {}

  SystemAssembler(TestSpaceType test, GridViewType grid_view_in)
    : BaseType(grid_view_in)
    , test_space_(test)
    , ansatz_space_(test)
  {}

  const TestSpaceType& test_space() const
  {
    return *test_space_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return *ansatz_space_;
  }

  using BaseType::add;

  template< class C, class M >
  void add(Spaces::ConstraintsInterface< C, RangeFieldType >& constraints,
           Stuff::LA::MatrixInterface< M, RangeFieldType >& matrix,
           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities< GridViewType >())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalMatrixConstraintsWrapper< TestSpaceType,
                                                     AnsatzSpaceType,
                                                     GridViewType,
                                                     typename C::derived_type,
                                                     typename M::derived_type > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_,
                                                        ansatz_space_,
                                                        where,
                                                        constraints.as_imp(),
                                                        matrix.as_imp()));
  } // ... add(...)

  /** \todo Investigate why the ConstraintsInterface is not used here! */
  template< class ConstraintsType, class V >
  void add(ConstraintsType& constraints,
           Stuff::LA::VectorInterface< V, RangeFieldType >& vector,
           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = vector.as_imp();
    assert(vector_imp.size() == test_space_->mapper().size());
    typedef internal::LocalVectorConstraintsWrapper< ThisType, ConstraintsType, VectorType > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, constraints, vector_imp));
  } // ... add(...)

  template< class L, class M >
  void add(const LocalAssembler::Codim0Matrix< L >& local_assembler,
           Stuff::LA::MatrixInterface< M, RangeFieldType >& matrix,
           const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = matrix.as_imp();
    assert(matrix_imp.rows() == test_space_->mapper().size());
    assert(matrix_imp.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalVolumeMatrixAssemblerWrapper< ThisType, LocalAssembler::Codim0Matrix< L >, MatrixType >
        WrapperType;
    this->codim0_functors_.emplace_back(
          new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  } // ... add(...)

  template< class Codim0Assembler, class M >
  void add_codim0_assembler(const Codim0Assembler& local_assembler,
                            Stuff::LA::MatrixInterface< M, RangeFieldType >& matrix,
                            const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = matrix.as_imp();
    assert(matrix_imp.rows() == test_space_->mapper().size());
    assert(matrix_imp.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalVolumeMatrixAssemblerWrapper< ThisType, Codim0Assembler, MatrixType > WrapperType;
    this->codim0_functors_.emplace_back(
          new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  } // ... add(...)

  template< class Codim0Assembler, class V >
  void add_codim0_assembler(const Codim0Assembler& local_assembler,
                            Stuff::LA::VectorInterface< V, RangeFieldType >& vector,
                            const ApplyOnWhichEntity* where = new DSG::ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = vector.as_imp();
    assert(vector_imp.size() == test_space_->mapper().size());
    typedef internal::LocalVolumeVectorAssemblerWrapper< ThisType, Codim0Assembler, VectorType > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  } // ... add(...)

  template< class L, class M >
  void add(const LocalAssembler::Codim1CouplingMatrix< L >& local_assembler,
           Stuff::LA::MatrixInterface< M, RangeFieldType >& matrix,
           const ApplyOnWhichIntersection* where
              = new DSG::ApplyOn::AllIntersections< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = matrix.as_imp();
    assert(matrix_imp.rows() == test_space_->mapper().size());
    assert(matrix_imp.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalFaceMatrixAssemblerWrapper< ThisType, LocalAssembler::Codim1CouplingMatrix< L >, MatrixType >
        WrapperType;
    this->codim1_functors_.emplace_back(
          new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  } // ... add(...)

  template< class L, class M >
  void add(const LocalAssembler::Codim1BoundaryMatrix< L >& local_assembler,
           Stuff::LA::MatrixInterface< M, RangeFieldType >& matrix,
           const ApplyOnWhichIntersection* where
              = new DSG::ApplyOn::AllIntersections< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = matrix.as_imp();
    assert(matrix_imp.rows() == test_space_->mapper().size());
    assert(matrix_imp.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalFaceMatrixAssemblerWrapper< ThisType, LocalAssembler::Codim1BoundaryMatrix< L >, MatrixType >
        WrapperType;
    this->codim1_functors_.emplace_back(
          new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  } // ... add(...)

  template< class L, class V >
  void add(const LocalAssembler::Codim0Vector< L >& local_assembler,
           Stuff::LA::VectorInterface< V, RangeFieldType >& vector,
           const ApplyOnWhichEntity* where
              = new DSG::ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = vector.as_imp();
    assert(vector_imp.size() == test_space_->mapper().size());
    typedef internal::LocalVolumeVectorAssemblerWrapper< ThisType, LocalAssembler::Codim0Vector< L >, VectorType >
        WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  } // ... add(...)

  template< class L, class V >
  void add(const LocalAssembler::Codim1Vector< L >& local_assembler,
           Stuff::LA::VectorInterface< V, RangeFieldType >& vector,
           const ApplyOnWhichIntersection* where
              = new DSG::ApplyOn::AllIntersections< GridViewType >())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast< VectorType& >(vector);
    assert(vector_imp.size() == test_space_->mapper().size());
    typedef internal::LocalFaceVectorAssemblerWrapper< ThisType, LocalAssembler::Codim1Vector< L >, VectorType >
        WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  } // ... add(...)

  void assemble(const bool use_tbb = false)
  {
    this->walk(use_tbb);
  }

private:
  const DS::PerThreadValue<const TestSpaceType> test_space_;
  const DS::PerThreadValue<const AnsatzSpaceType> ansatz_space_;
}; // class SystemAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
