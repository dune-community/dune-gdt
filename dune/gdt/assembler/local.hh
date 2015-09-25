// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMLBER_LOCAL_HH
#define DUNE_GDT_ASSEMLBER_LOCAL_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/grid/walker/wrapper.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localoperator/interfaces.hh>
#include <dune/gdt/localfunctional/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


template< class LocalVolumeTwoFormType >
class LocalVolumeTwoFormAssembler
{
  static_assert(is_local_volume_twoform< LocalVolumeTwoFormType >::value,
                "LocalVolumeTwoFormType has to be derived from LocalVolumeTwoFormInterface!");

public:
  explicit LocalVolumeTwoFormAssembler(const LocalVolumeTwoFormType& local_twoform)
    : local_volume_twoform_(local_twoform)
  {}

  /**
   *  \tparam T           Traits of the SpaceInterface implementation, representing the type of test_space
   *  \tparam A           Traits of the SpaceInterface implementation, representing the type of ansatz_space
   *  \tparam *d          dimDomain of test_space (* == T) or ansatz_space (* == A)
   *  \tparam *r          dimRange of test_space (* == T) or ansatz_space (* == A)
   *  \tparam *rC         dimRangeCols of test_space (* == T) or ansatz_space (* == A)
   *  \tparam EntityType  A model of Dune::Entity< 0 >
   *  \tparam M           Traits of the Dune::Stuff::LA::Container::MatrixInterface implementation, representing the type of global_matrix
   *  \tparam R           RangeFieldType, i.e. double
   */
  template< class T, size_t Td, size_t Tr, size_t TrC, class A, size_t Ad, size_t Ar, size_t ArC, class EntityType, class M, class R >
  void assemble(const SpaceInterface< T, Td, Tr, TrC >& test_space,
                const SpaceInterface< A, Ad, Ar, ArC >& ansatz_space,
                const EntityType& entity,
                Stuff::LA::MatrixInterface< M, R >& global_matrix) const
  {
    // prepare
    const size_t rows = test_space.mapper().numDofs(entity);
    const size_t cols = ansatz_space.mapper().numDofs(entity);
    Dune::DynamicMatrix< R > local_matrix(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // apply local two-form
    const auto test_base   = test_space.base_function_set(entity);
    const auto ansatz_base = ansatz_space.base_function_set(entity);
    assert(test_base.size()   == rows);
    assert(ansatz_base.size() == cols);
    local_volume_twoform_.apply2(test_base, ansatz_base, local_matrix);
    // write local matrix to global
    const auto global_row_indices = test_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    const auto global_col_indices = ansatz_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    assert(global_row_indices.size() == rows);
    assert(global_col_indices.size() == cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& local_matrix_row = local_matrix[ii];
      const size_t global_ii = global_row_indices[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const size_t global_jj = global_col_indices[jj];
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix_row[jj]);
      }
    } // write local matrix to global
  } // ... assemble(...)

private:
  const LocalVolumeTwoFormType& local_volume_twoform_;
}; // class LocalVolumeTwoFormAssembler


template< class GridViewImp, class LocalVolumeTwoFormType, class TestFunctionType, class AnsatzFunctionType, class FieldType >
class LocalVolumeTwoFormAccumulator
  : public Stuff::Grid::internal::Codim0ReturnObject< GridViewImp, FieldType >
{
  static_assert(std::is_base_of< LocalVolumeTwoFormInterface< typename LocalVolumeTwoFormType::Traits >,
                                 LocalVolumeTwoFormType >::value,
                "LocalVolumeTwoFormType has to be derived from LocalVolumeTwoFormInterface!");
  static_assert(Stuff::is_localizable_function< TestFunctionType >::value,
                "TestFunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function< AnsatzFunctionType >::value,
                "AnsatzFunctionType has to be derived from Stuff::LocalizableFunctionInterface!");

  typedef LocalVolumeTwoFormAccumulator
      < GridViewImp, LocalVolumeTwoFormType, TestFunctionType, AnsatzFunctionType, FieldType > ThisType;
  typedef Stuff::Grid::internal::Codim0ReturnObject< GridViewImp, FieldType >                  BaseType;
public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType   EntityType;

  LocalVolumeTwoFormAccumulator(const GridViewType& grd_vw,
                                const LocalVolumeTwoFormType& local_op,
                                const TestFunctionType& test_function,
                                const AnsatzFunctionType& ansatz_function,
                                const Stuff::Grid::ApplyOn::WhichEntity< GridViewType >& where)
    : grid_view_(grd_vw)
    , local_operator_(local_op)
    , test_function_(test_function)
    , ansatz_function_(ansatz_function)
    , result_(0)
    , finalized_(false)
    , where_(where)
  {}

  LocalVolumeTwoFormAccumulator(const ThisType& other) = default;
  virtual ~LocalVolumeTwoFormAccumulator()             = default;

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const override final
  {
    return where_.apply_on(grid_view, entity);
  }

  virtual FieldType compute_locally(const EntityType& entity) override final
  {
    DynamicMatrix< FieldType > local_twoform_result(1, 1, 0.); // \todo: make mutable member, after SMP refactor
    this->local_operator_.apply2(*test_function_.local_function(entity),
                                 *ansatz_function_.local_function(entity),
                                 local_twoform_result);
    return local_twoform_result[0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType& entity) override
  {
    *result_ += compute_locally(entity);
  }

  virtual void finalize() override
  {
    if (!finalized_) {
      finalized_result_ = result_.sum();
      finalized_result_ = grid_view_.comm().sum(finalized_result_);
      finalized_ = true;
    }
  } // ... finalize(...)

  virtual FieldType result() const override final
  {
    if (!finalized_)
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong, "Call finalize() first!");
    return finalized_result_;
  }

private:
  const GridViewType& grid_view_;
  const LocalVolumeTwoFormType& local_operator_;
  const TestFunctionType& test_function_;
  const AnsatzFunctionType& ansatz_function_;
  DS::PerThreadValue< FieldType > result_;
  bool finalized_;
  const Stuff::Grid::ApplyOn::WhichEntity< GridViewType >& where_;
  FieldType finalized_result_;
}; // class LocalVolumeTwoFormAccumulator


template< class GridViewType, class LocalOperatorType, class SourceType, class RangeType >
class LocalOperatorApplicator
  : public Stuff::Grid::internal::Codim0Object< GridViewType >
{
  static_assert(is_local_operator< LocalOperatorType >::value,
                "LocalOperatorType has to be derived from LocalOperatorInterface!");
  static_assert(Stuff::is_localizable_function< SourceType >::value,
                "SourceType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_discrete_function< RangeType >::value,
                "RangeType has to be a DiscreteFunctionDefault!");
  typedef Stuff::Grid::internal::Codim0Object< GridViewType > BaseType;
public:
  using typename BaseType::EntityType;

  LocalOperatorApplicator(const GridViewType& grid_view,
                          const LocalOperatorType& local_operator,
                          const SourceType& source,
                          RangeType& range,
                          const Stuff::Grid::ApplyOn::WhichEntity< GridViewType >& where)
    : grid_view_(grid_view)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {}

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const
  {
    return where_.apply_on(grid_view, entity);
  }

  virtual void apply_local(const EntityType& entity)
  {
    local_operator_.apply(source_, *range_.local_discrete_function(entity));
  }

private:
  const GridViewType& grid_view_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const Stuff::Grid::ApplyOn::WhichEntity< GridViewType >& where_;
}; // class LocalOperatorApplicator


template< class LocalVolumeFunctionalType >
class LocalVolumeFunctionalAssembler
{
  static_assert(is_local_volume_functional< LocalVolumeFunctionalType >::value,
                "LocalVolumeFunctionalType has to be derived from LocalVolumeFunctionalInterface!");

public:
  explicit LocalVolumeFunctionalAssembler(const LocalVolumeFunctionalType& local_volume_functional)
    : local_volume_functional_(local_volume_functional)
  {}

  /**
   *  \tparam S          Traits of the SpaceInterface implementation, representing the type of test_space
   *  \tparam d          dimDomain of test_space
   *  \tparam r          dimRange of test_space
   *  \tparam rC         dimRangeCols of test_space
   *  \tparam EntityType A model of Dune::Entity< 0 >
   *  \tparam V          Traits of the Dune::Stuff::LA::Container::VectorInterface implementation, representing the type of global_vector
   *  \tparam R          RangeFieldType, i.e. double
   */
  template< class S, size_t d, size_t r, size_t rC, class EntityType, class V, class R >
  void assemble(const SpaceInterface< S, d, r, rC >& test_space,
                const EntityType& entity,
                Stuff::LA::VectorInterface< V, R >& global_vector) const
  {
    // prepare
    const size_t size = test_space.mapper().numDofs(entity);
    Dune::DynamicVector< R > local_vector(size, 0.); // \todo: make mutable member, after SMP refactor
    // apply local functional
    const auto test_basis = test_space.base_function_set(entity);
    assert(test_basis.size() == size);
    local_volume_functional_.apply(test_basis, local_vector);
    // write local vector to global
    const auto global_indices = test_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    assert(global_indices.size() == size);
    for (size_t jj = 0; jj < size; ++jj)
      global_vector.add_to_entry(global_indices[jj], local_vector[jj]);
  } // ... assemble(...)

private:
  const LocalVolumeFunctionalType& local_volume_functional_;
}; // class LocalVolumeFunctionalAssembler


template< class LocalFunctionalType >
class LocalFaceFunctionalAssembler
{
  static_assert(std::is_base_of< LocalFaceFunctionalInterface< typename LocalFunctionalType::Traits >,
                                 LocalFunctionalType >::value,
                "LocalFunctionalType has to be derived from LocalFaceFunctionalInterface!");
public:
  explicit LocalFaceFunctionalAssembler(const LocalFunctionalType& local_face_functional)
    : local_face_functional_(local_face_functional)
  {}

  template< class T, size_t d, size_t r, size_t rC, class IntersectionType, class V, class R >
  void assemble(const SpaceInterface< T, d, r, rC >& test_space,
                     const IntersectionType& intersection,
                     Stuff::LA::VectorInterface< V, R >& global_vector) const
  {
    // prepare
    const auto entity_ptr = intersection.inside();
    const auto& entity = *entity_ptr;
    const size_t size = test_space.mapper().numDofs(entity);
    Dune::DynamicVector< R > local_vector(size, 0.); // \todo: make mutable member, after SMP refactor
    // apply local functional
    const auto test_basis = test_space.base_function_set(entity);
    assert(test_basis.size() == size);
    local_face_functional_.apply(test_basis, intersection, local_vector);
    // write local vector to global
    const auto global_indices = test_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    assert(global_indices.size() == size);
    for (size_t jj = 0; jj < size; ++jj)
      global_vector.add_to_entry(global_indices[jj], local_vector[jj]);
  } // ... assemble(...)

private:
  const LocalFunctionalType& local_face_functional_;
}; // class LocalFaceFunctionalAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMLBER_LOCAL_HH
