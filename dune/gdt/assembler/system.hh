// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_HH

#include <type_traits>
#include <memory>

#if DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) //&& HAVE_TBB //EXADUNE
#include <dune/grid/utility/partitioning/seedlist.hh>
#endif

#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/common/parallel/partitioner.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>

#include "local/codim0.hh"
#include "local/codim1.hh"
#include "wrapper.hh"

namespace Dune {
namespace GDT {


template <class TestSpaceImp, class GridViewImp = typename TestSpaceImp::GridViewType,
          class AnsatzSpaceImp                  = TestSpaceImp>
class SystemAssembler : public Stuff::Grid::Walker<GridViewImp>
{
  static_assert(std::is_base_of<SpaceInterface<typename TestSpaceImp::Traits>, TestSpaceImp>::value,
                "TestSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename AnsatzSpaceImp::Traits>, AnsatzSpaceImp>::value,
                "AnsatzSpaceImp has to be derived from SpaceInterface!");
  typedef Stuff::Grid::Walker<GridViewImp> BaseType;
  typedef SystemAssembler<TestSpaceImp, GridViewImp, AnsatzSpaceImp> ThisType;

public:
  typedef TestSpaceImp TestSpaceType;
  typedef AnsatzSpaceImp AnsatzSpaceType;

  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  typedef Stuff::Grid::ApplyOn::WhichEntity<GridViewType> ApplyOnWhichEntity;
  typedef Stuff::Grid::ApplyOn::WhichIntersection<GridViewType> ApplyOnWhichIntersection;

  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz, const GridViewType& grid_view)
    : BaseType(grid_view)
    , test_space_(test)
    , ansatz_space_(ansatz)
  {
  }

  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz)
    : BaseType(*(test.grid_view()))
    , test_space_(test)
    , ansatz_space_(ansatz)
  {
  }

  SystemAssembler(const TestSpaceType& test)
    : BaseType(*(test.grid_view()))
    , test_space_(test)
    , ansatz_space_(test_space_)
  {
  }

  SystemAssembler(const TestSpaceType& test, const GridViewType& grid_view)
    : BaseType(grid_view)
    , test_space_(test)
    , ansatz_space_(test_space_)
  {
  }

  const TestSpaceType& test_space() const
  {
    return test_space_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return ansatz_space_;
  }

  using BaseType::add;

  template <class ConstraintsType, class M>
  void add(ConstraintsType& constraints, Stuff::LA::MatrixInterface<M>& matrix,
           const ApplyOnWhichEntity* where = new Stuff::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast<MatrixType&>(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef internal::LocalMatrixConstraintsWrapper<ThisType, ConstraintsType, MatrixType> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, constraints, matrix_imp));
  } // ... add(...)

  template <class ConstraintsType, class V>
  void add(ConstraintsType& constraints, Stuff::LA::VectorInterface<V>& vector,
           const ApplyOnWhichEntity* where = new Stuff::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast<VectorType&>(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef internal::LocalVectorConstraintsWrapper<ThisType, ConstraintsType, VectorType> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, constraints, vector_imp));
  } // ... add(...)

  template <class L, class M>
  void add(const LocalAssembler::Codim0Matrix<L>& local_assembler, Stuff::LA::MatrixInterface<M>& matrix,
           const ApplyOnWhichEntity* where = new Stuff::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast<MatrixType&>(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef internal::LocalVolumeMatrixAssemblerWrapper<ThisType, LocalAssembler::Codim0Matrix<L>, MatrixType>
        WrapperType;
    this->codim0_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  } // ... add(...)

  template <class Codim0Assembler, class M>
  void add_codim0_assembler(const Codim0Assembler& local_assembler, Stuff::LA::MatrixInterface<M>& matrix,
                            const ApplyOnWhichEntity* where = new Stuff::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast<MatrixType&>(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef internal::LocalVolumeMatrixAssemblerWrapper<ThisType, Codim0Assembler, MatrixType> WrapperType;
    this->codim0_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  } // ... add(...)

  template <class Codim0Assembler, class V>
  void add_codim0_assembler(const Codim0Assembler& local_assembler, Stuff::LA::VectorInterface<V>& vector,
                            const ApplyOnWhichEntity* where = new Stuff::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast<VectorType&>(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef internal::LocalVolumeVectorAssemblerWrapper<ThisType, Codim0Assembler, VectorType> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  } // ... add(...)

  template <class L, class M>
  void add(const LocalAssembler::Codim1CouplingMatrix<L>& local_assembler, Stuff::LA::MatrixInterface<M>& matrix,
           const ApplyOnWhichIntersection* where = new Stuff::Grid::ApplyOn::AllIntersections<GridViewType>())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast<MatrixType&>(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef internal::LocalFaceMatrixAssemblerWrapper<ThisType, LocalAssembler::Codim1CouplingMatrix<L>, MatrixType>
        WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  } // ... add(...)

  template <class L, class M>
  void add(const LocalAssembler::Codim1BoundaryMatrix<L>& local_assembler, Stuff::LA::MatrixInterface<M>& matrix,
           const ApplyOnWhichIntersection* where = new Stuff::Grid::ApplyOn::AllIntersections<GridViewType>())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast<MatrixType&>(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef internal::LocalFaceMatrixAssemblerWrapper<ThisType, LocalAssembler::Codim1BoundaryMatrix<L>, MatrixType>
        WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  } // ... add(...)

  template <class L, class V>
  void add(const LocalAssembler::Codim0Vector<L>& local_assembler, Stuff::LA::VectorInterface<V>& vector,
           const ApplyOnWhichEntity* where = new Stuff::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast<VectorType&>(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef internal::LocalVolumeVectorAssemblerWrapper<ThisType, LocalAssembler::Codim0Vector<L>, VectorType>
        WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  } // ... add(...)

  template <class L, class V>
  void add(const LocalAssembler::Codim1Vector<L>& local_assembler, Stuff::LA::VectorInterface<V>& vector,
           const ApplyOnWhichIntersection* where = new Stuff::Grid::ApplyOn::AllIntersections<GridViewType>())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast<VectorType&>(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef internal::LocalFaceVectorAssemblerWrapper<ThisType, LocalAssembler::Codim1Vector<L>, VectorType>
        WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  } // ... add(...)

  void assemble(const bool clear_stack = true)
  {
    this->walk(clear_stack);
  }

#if DUNE_VERSION_NEWER(DUNE_COMMON, 3, 9) //&& HAVE_TBB //EXADUNE
  void tbb_assemble(const bool clear_stack = true)
  {
    Stuff::IndexSetPartitioner<GridViewType> partioner(this->grid_view_.indexSet());
    SeedListPartitioning<typename GridViewType::Grid, 0> partitioning(this->grid_view_, partioner);
    this->tbb_walk(partitioning, clear_stack);
  }
#endif

private:
  const TestSpaceType& test_space_;
  const AnsatzSpaceType& ansatz_space_;
}; // class SystemAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
