// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   René Fritze     (2016 - 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2016 - 2018)

#ifndef DUNE_GDT_OPERATORS_MATRIX_BASED_OPERATOR_HH
#define DUNE_GDT_OPERATORS_MATRIX_BASED_OPERATOR_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/matrix-interface.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/assembler/bilinear-form-assemblers.hh>
#include <dune/gdt/local/assembler/operator-fd-jacobian-assemblers.hh>
#include <dune/gdt/local/bilinear-forms/interfaces.hh>
#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Base class for linear operators which are given by an assembled matrix.
 *
 * \note See OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 * \sa MatrixOperator
 */
template <class M, class SGV, size_t s_r = 1, size_t s_rC = 1, size_t r_r = s_r, size_t r_rC = s_rC, class RGV = SGV>
class ConstMatrixOperator : public OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>
{
  using ThisType = ConstMatrixOperator;
  using BaseType = OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;

public:
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  ConstMatrixOperator(const SourceSpaceType& source_spc, const RangeSpaceType& range_spc, const MatrixType& mat)
    : source_space_(source_spc)
    , range_space_(range_spc)
    , matrix_(mat)
    , linear_solver_(matrix_, source_space_.dof_communicator())
  {
    DUNE_THROW_IF(matrix_.rows() != range_space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "matrix_.rows() = " << matrix_.rows()
                                      << "\n   range_space_.mapper().size() = " << range_space_.mapper().size());
    DUNE_THROW_IF(matrix_.cols() != source_space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "matrix_.cols() = " << matrix_.cols()
                                      << "\n   source_space_.mapper().size() = " << source_space_.mapper().size());
  } // ConstMatrixOperator(...)

  ConstMatrixOperator(ThisType&& source)
    : source_space_(source.source_space_)
    , range_space_(source.range_space_)
    , matrix_(source.matrix_)
    , linear_solver_(matrix_, source_space_.dof_communicator())
  {}

  bool linear() const override final
  {
    return true;
  }

  const SourceSpaceType& source_space() const override
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override
  {
    return range_space_;
  }

  const MatrixType& matrix() const
  {
    return matrix_;
  }

  using BaseType::apply;

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& /*param*/ = {}) const override
  {
    try {
      matrix_.mv(source, range);
    } catch (const XT::Common::Exceptions::shapes_do_not_match& ee) {
      DUNE_THROW(Exceptions::operator_error,
                 "when applying matrix to source and range dofs!\n\nThis was the original error: " << ee.what());
    }
  } // ... apply(...)

  std::vector<std::string> invert_options() const override
  {
    return linear_solver_.types();
  }

  XT::Common::Configuration invert_options(const std::string& type) const override
  {
    try {
      return linear_solver_.options(type);
    } catch (const XT::Common::Exceptions::configuration_error& ee) {
      DUNE_THROW(Exceptions::operator_error,
                 "when accessing linear solver options!\n\nThis was the original error: " << ee.what());
    }
  } // ... invert_options(...)

  using BaseType::apply_inverse;

  void apply_inverse(const VectorType& range,
                     VectorType& source,
                     const XT::Common::Configuration& opts,
                     const XT::Common::Parameter& /*param*/ = {}) const override
  {
    try {
      linear_solver_.apply(range, source, opts);
    } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
      DUNE_THROW(Exceptions::operator_error,
                 "when applying linear solver!\n\nThis was the original error: " << ee.what());
    }
  } // ... apply_inverse(...)

  std::vector<std::string> jacobian_options() const override
  {
    return {"self"};
  }

  XT::Common::Configuration jacobian_options(const std::string& type) const override
  {
    DUNE_THROW_IF(type != jacobian_options().at(0),
                  Exceptions::operator_error,
                  "requested jacobian type is not one of the available ones!\n\n"
                      << "type = " << type << "\njacobian_options() = " << jacobian_options());
    return {{"type", jacobian_options().at(0)}};
  } // ... jacobian_options(...)

  using BaseType::jacobian;

  void jacobian(const VectorType& /*source*/,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& /*param*/ = {}) const override
  {
    DUNE_THROW_IF(
        !opts.has_key("type"), Exceptions::operator_error, "missing key 'type' in given opts!\n\nopts = " << opts);
    const auto type = opts.get<std::string>("type");
    DUNE_THROW_IF(type != jacobian_options().at(0),
                  Exceptions::operator_error,
                  "requested jacobian type is not one of the available ones!\n\n"
                      << "type = " << type << "\njacobian_options() = " << jacobian_options());
    jacobian_op.matrix() += matrix_;
  } // ... jacobian(...)

private:
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const MatrixType& matrix_;
  const XT::LA::Solver<MatrixType, typename SourceSpaceType::DofCommunicatorType> linear_solver_;
}; // class ConstMatrixOperator


template <class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC, class M>
ConstMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, SGV, s_r, s_rC, r_r, r_rC, RGV>
make_matrix_operator(const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                     const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                     const XT::LA::MatrixInterface<M>& matrix)
{
  return ConstMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, SGV, s_r, s_rC, r_r, r_rC, RGV>(
      source_space, range_space, matrix.as_imp());
}


template <class GV, size_t r, size_t rC, class F, class M>
ConstMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC>
make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space, const XT::LA::MatrixInterface<M>& matrix)
{
  return ConstMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC>(
      space, space, matrix.as_imp());
}


/**
 * \brief Base class for linear operators which are assembled into a matrix.
 *
 * We derive from the XT::Grid::Walker and povide custom append() methods to allow to add local element and intersection
 * operators. In contrast to the GlobalAssembler we already hold the target matrix (or create one of appropriate size
 * and given sparsity pattern), into which we want to assemble. The operator is assembled by walking over the given
 * assembly_grid_view, if you want to assemble an operator only on a smaller grid view, consider to append this operator
 * to another walker or provide appropriate filters.
 *
 * \note See ConstMatrixOperator and OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 * \sa ConstMatrixOperator
 * \sa XT::Grid::Walker
 */
template <class M, class SGV, size_t s_r = 1, size_t s_rC = 1, size_t r_r = s_r, size_t r_rC = s_rC, class RGV = SGV>
class MatrixOperator
  : XT::Common::StorageProvider<M>
  , public ConstMatrixOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>
  , public XT::Grid::Walker<SGV>
{
  using ThisType = MatrixOperator;
  using MatrixStorage = XT::Common::StorageProvider<M>;
  using OperatorBaseType = ConstMatrixOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;
  using WalkerBaseType = XT::Grid::Walker<SGV>;

public:
  using AssemblyGridViewType = SGV;

  using typename OperatorBaseType::F;
  using typename OperatorBaseType::FieldType;
  using typename OperatorBaseType::MatrixOperatorType;
  using typename OperatorBaseType::MatrixType;
  using typename OperatorBaseType::RangeSpaceType;
  using typename OperatorBaseType::SourceSpaceType;
  using typename OperatorBaseType::V;
  using typename OperatorBaseType::VectorType;

  using typename WalkerBaseType::ElementType;
  using ElementFilterType = XT::Grid::ElementFilter<SGV>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<SGV>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<SGV>;
  using ApplyOnAllIntersections = XT::Grid::ApplyOn::AllIntersections<SGV>;

  using typename WalkerBaseType::E;
  using typename WalkerBaseType::I;

  /**
   * Ctor which accept an existing matrix into which to assemble.
   */
  MatrixOperator(AssemblyGridViewType assembly_grid_view,
                 const SourceSpaceType& source_spc,
                 const RangeSpaceType& range_spc,
                 MatrixType& mat)
    : MatrixStorage(mat)
    , OperatorBaseType(source_spc, range_spc, MatrixStorage::access())
    , WalkerBaseType(assembly_grid_view)
    , scaling(1.)
    , assembled_(false)
  {
    // to detect assembly
    this->append(
        [](/*prepare nothing*/) {}, [](const auto&) { /*apply nothing*/ }, [&](/*finalize*/) { assembled_ = true; });
  }

  /**
   * Ctor which creates an appropriate matrix into which to assemble from a given sparsity pattern.
   */

  MatrixOperator(AssemblyGridViewType assembly_grid_view,
                 const SourceSpaceType& source_spc,
                 const RangeSpaceType& range_spc,
                 const XT::LA::SparsityPatternDefault& pattern)
    : MatrixStorage(new MatrixType(range_spc.mapper().size(), source_spc.mapper().size(), pattern))
    , OperatorBaseType(source_spc, range_spc, MatrixStorage::access())
    , WalkerBaseType(assembly_grid_view)
    , scaling(1.)
    , assembled_(false)
  {
    // to detect assembly
    this->append(
        [](/*prepare nothing*/) {}, [](const auto&) { /*apply nothing*/ }, [&](/*finalize*/) { assembled_ = true; });
  }

  using OperatorBaseType::matrix;

  MatrixType& matrix()
  {
    return MatrixStorage::access();
  }

  FieldType scaling;

  void clear()
  {
    WalkerBaseType::clear();
    assembled_ = false;
  }

  using WalkerBaseType::append;

  ThisType& append(const LocalElementBilinearFormInterface<E, r_r, r_rC, F, F, s_r, s_rC, F>& local_bilinear_form,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    using LocalAssemblerType = LocalElementBilinearFormAssembler<M, SGV, r_r, r_rC, F, RGV, SGV, s_r, s_rC, F>;
    this->append(new LocalAssemblerType(this->range_space(),
                                        this->source_space(),
                                        local_bilinear_form,
                                        MatrixStorage::access(),
                                        param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
                 filter);
    return *this;
  }

  ThisType& append(const LocalIntersectionBilinearFormInterface<I, r_r, r_rC, F, F, s_r, s_rC, F>& local_bilinear_form,
                   const XT::Common::Parameter& param = {},
                   const IntersectionFilterType& filter = ApplyOnAllIntersections())
  {
    using LocalAssemblerType =
        LocalIntersectionBilinearFormAssembler<MatrixType, AssemblyGridViewType, r_r, r_rC, F, RGV, SGV, s_r, s_rC, F>;
    this->append(new LocalAssemblerType(this->range_space(),
                                        this->source_space(),
                                        local_bilinear_form,
                                        MatrixStorage::access(),
                                        param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
                 filter);
    return *this;
  } // ... append(...)

  /// \{
  /// \name Variants to compute the jacobian of the appended local operator by finite differences.

  ThisType& append(const LocalElementOperatorInterface<V, SGV, s_r, s_rC, F, r_r, r_rC>& local_operator,
                   const VectorType& source,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    this->append(new LocalElementOperatorFiniteDifferenceJacobianAssembler<M, SGV, s_r, s_rC, F, r_r, r_rC>(
                     this->source_space(),
                     this->range_space(),
                     MatrixStorage::access(),
                     source,
                     local_operator,
                     param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
                 filter);
    return *this;
  }

  ThisType& append(const LocalIntersectionOperatorInterface<I, V, SGV, s_r, s_rC, F, r_r, r_rC>& local_operator,
                   const VectorType& source,
                   const XT::Common::Parameter& param = {},
                   const IntersectionFilterType& filter = ApplyOnAllIntersections())
  {
    this->append(new LocalIntersectionOperatorFiniteDifferenceJacobianAssembler<M, SGV, s_r, s_rC, F, r_r, r_rC>(
                     this->source_space(),
                     this->range_space(),
                     MatrixStorage::access(),
                     source,
                     local_operator,
                     param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
                 filter);
    return *this;
  }

  /// \}

  ThisType& assemble(const bool use_tbb = false) override final
  {
    if (!assembled_) {
      // This clears all appended operators, which is ok, since we are done after assembling once!
      this->walk(use_tbb);
      assembled_ = true;
    }
    return *this;
  } // ... assemble(...)

  template <class EntityRange>
  ThisType& assemble_range(const EntityRange& entity_range)
  {
    if (!assembled_) {
      // This clears all appended operators, which is ok, since we are done after assembling once!
      this->walk_range(entity_range);
      assembled_ = true;
    }
    return *this;
  } // ... assemble(...)

  using OperatorBaseType::jacobian;

  /// \todo Store appended local bilinear forms and append them to jacobian_op?
  void jacobian(const VectorType& source,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& param = {}) const override
  {
    DUNE_THROW_IF(!assembled_, Exceptions::operator_error, "This operator has to be assembled to provide a jacobian!");
    OperatorBaseType::jacobian(source, jacobian_op, opts, param);
  }

private:
  bool assembled_;
}; // class MatrixOperator


/// \name Variants of make_matrix_operator for a given matrix.
/// \{

template <class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC, class M>
MatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, SGV, s_r, s_rC, r_r, r_rC, RGV>
make_matrix_operator(SGV assembly_grid_view,
                     const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                     const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                     XT::LA::MatrixInterface<M>& matrix)
{
  return MatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, SGV, s_r, s_rC, r_r, r_rC, RGV>(
      assembly_grid_view, source_space, range_space, matrix.as_imp());
}

template <class GV, size_t r, size_t rC, class F, class M>
MatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC> make_matrix_operator(
    GV assembly_grid_view, const SpaceInterface<GV, r, rC, F>& space, XT::LA::MatrixInterface<M>& matrix)
{
  return MatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC>(
      assembly_grid_view, space, space, matrix.as_imp());
}

template <class GV, size_t r, size_t rC, class F, class M>
MatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC>
make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space, XT::LA::MatrixInterface<M>& matrix)
{
  return MatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC>(
      space.grid_view(), space, space, matrix.as_imp());
}

/// \}
/// \name Variants of make_matrix_operator, where an appropriate matrix is created from a given pattern
/// \{

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(assembly_grid_view, source_space, range_space, pattern);
\endcode
 */
template <class MatrixType, class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
                        MatrixOperator<MatrixType, SGV, s_r, s_rC, r_r, r_rC, RGV>>::type
make_matrix_operator(SGV assembly_grid_view,
                     const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                     const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                     const XT::LA::SparsityPatternDefault& pattern)
{
  return MatrixOperator<MatrixType, SGV, s_r, s_rC, r_r, r_rC, RGV>(
      assembly_grid_view, source_space, range_space, pattern);
}

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(assembly_grid_view, space, pattern);
\endcode
 */
template <class MatrixType, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value, MatrixOperator<MatrixType, GV, r, rC>>::type
make_matrix_operator(GV assembly_grid_view,
                     const SpaceInterface<GV, r, rC, F>& space,
                     const XT::LA::SparsityPatternDefault& pattern)
{
  return MatrixOperator<MatrixType, GV, r, rC>(assembly_grid_view, space, space, pattern);
}

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(space, pattern);
\endcode
 */
template <class MatrixType, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value, MatrixOperator<MatrixType, GV, r, rC>>::type
make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space, const XT::LA::SparsityPatternDefault& pattern)
{
  return MatrixOperator<MatrixType, GV, r, rC>(space.grid_view(), space, space, pattern);
}

/// \}
/// \name Variants of make_matrix_operator, where an appropriate matrix is created from given stencil
/// \{

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(source_space, range_space, stencil);
\endcode
 */
template <class MatrixType, class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
                        MatrixOperator<MatrixType, SGV, s_r, s_rC, r_r, r_rC, RGV>>::type
make_matrix_operator(SGV assembly_grid_view,
                     const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                     const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                     const Stencil stencil = Stencil::element_and_intersection)
{
  return MatrixOperator<MatrixType, SGV, s_r, s_rC, r_r, r_rC, RGV>(
      assembly_grid_view,
      source_space,
      range_space,
      make_sparsity_pattern(range_space, source_space, assembly_grid_view, stencil));
}

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(assembly_grid_view, space, stencil);
\endcode
 */
template <class MatrixType, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value, MatrixOperator<MatrixType, GV, r, rC>>::type
make_matrix_operator(GV assembly_grid_view,
                     const SpaceInterface<GV, r, rC, F>& space,
                     const Stencil stencil = Stencil::element_and_intersection)
{
  return MatrixOperator<MatrixType, GV, r, rC>(
      assembly_grid_view, space, space, make_sparsity_pattern(space, assembly_grid_view, stencil));
}

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(space, stencil);
\endcode
 */
template <class MatrixType, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value, MatrixOperator<MatrixType, GV, r, rC>>::type
make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space,
                     const Stencil stencil = Stencil::element_and_intersection)
{
  return MatrixOperator<MatrixType, GV, r, rC>(space.grid_view(), space, space, make_sparsity_pattern(space, stencil));
}

/// \}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_MATRIX_BASED_OPERATOR_HH
