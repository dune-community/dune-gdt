// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

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
#include <dune/gdt/local/assembler/two-form-assemblers.hh>
#include <dune/gdt/local/bilinear-forms/interfaces.hh>
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
 * \sa MatrixBasedOperator
 */
template <class Matrix,
          class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          class F = double,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RF = double,
          class RGV = SGV,
          class SourceVector = typename XT::LA::Container<typename Matrix::ScalarType, Matrix::vector_type>::VectorType,
          class RangeVector = SourceVector>
class ConstMatrixBasedOperator
    : public OperatorInterface<SourceVector, SGV, s_r, s_rC, SF, F, r_r, r_rC, RF, RGV, RangeVector>
{
  static_assert(XT::LA::is_matrix<Matrix>::value, "");

  using ThisType =
      ConstMatrixBasedOperator<Matrix, SGV, s_r, s_rC, SF, F, r_r, r_rC, RF, RGV, SourceVector, RangeVector>;
  using BaseType = OperatorInterface<SourceVector, SGV, s_r, s_rC, SF, F, r_r, r_rC, RF, RGV, RangeVector>;

public:
  using MatrixType = Matrix;
  using DofFieldType = typename MatrixType::ScalarType;

  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;

  using typename BaseType::SourceVectorType;
  using typename BaseType::RangeVectorType;

  ConstMatrixBasedOperator(const SourceSpaceType& source_spc, const RangeSpaceType& range_spc, const MatrixType& mat)
    : source_space_(source_spc)
    , range_space_(range_spc)
    , matrix_(mat)
    , linear_solver_(matrix_, source_space_.dof_communicator())
  {
    DUNE_THROW_IF(matrix_.rows() != range_space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "matrix_.rows() = " << matrix_.rows() << "\n   range_space_.mapper().size() = "
                                      << range_space_.mapper().size());
    DUNE_THROW_IF(matrix_.cols() != source_space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "matrix_.cols() = " << matrix_.cols() << "\n   source_space_.mapper().size() = "
                                      << source_space_.mapper().size());
  } // ConstMatrixBasedOperator(...)

  ConstMatrixBasedOperator(ThisType&& source)
    : source_space_(source.source_space_)
    , range_space_(source.range_space_)
    , matrix_(source.matrix_)
    , linear_solver_(matrix_, source_space_.dof_communicator())
  {
  }

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

  void apply(const SourceVectorType& source,
             RangeVectorType& range,
             const XT::Common::Parameter& /*param*/ = {}) const override
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

  void apply_inverse(const RangeVectorType& range,
                     SourceVectorType& source,
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
    return {"direct"};
  }

  XT::Common::Configuration jacobian_options(const std::string& type) const override
  {
    DUNE_THROW_IF(type != jacobian_options().at(0),
                  Exceptions::operator_error,
                  "requested jacobian type is not one of the available ones!\n\n"
                      << "type = "
                      << type
                      << "\njacobian_options() = "
                      << jacobian_options());
    return {{"type", jacobian_options().at(0)}};
  } // ... jacobian_options(...)

  using BaseType::jacobian;

  std::shared_ptr<BaseType> jacobian(const SourceVectorType& /*source*/,
                                     const XT::Common::Configuration& opts,
                                     const XT::Common::Parameter& /*param*/ = {}) const override
  {
    DUNE_THROW_IF(
        !opts.has_key("type"), Exceptions::operator_error, "missing key 'type' in given opts!\n\nopts = " << opts);
    const auto type = opts.get<std::string>("type");
    DUNE_THROW_IF(type != jacobian_options().at(0),
                  Exceptions::operator_error,
                  "requested jacobian type is not one of the available ones!\n\n"
                      << "type = "
                      << type
                      << "\njacobian_options() = "
                      << jacobian_options());
    return std::make_shared<ThisType>(source_space_, range_space_, matrix_);
  } // ... jacobian(...)

private:
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const MatrixType& matrix_;
  const XT::LA::Solver<MatrixType, typename SourceSpaceType::DofCommunicatorType> linear_solver_;
}; // class ConstMatrixBasedOperator


template <class SGV, size_t s_r, size_t s_rC, class SF, class RGV, size_t r_r, size_t r_rC, class RF, class M>
ConstMatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type,
                         SGV,
                         s_r,
                         s_rC,
                         SF,
                         typename XT::Common::multiplication_promotion<SF, RF>::type,
                         r_r,
                         r_rC,
                         RF,
                         RGV>
make_matrix_operator(const SpaceInterface<SGV, s_r, s_rC, SF>& source_space,
                     const SpaceInterface<RGV, r_r, r_rC, RF>& range_space,
                     const XT::LA::MatrixInterface<M>& matrix)
{
  return ConstMatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type,
                                  SGV,
                                  s_r,
                                  s_rC,
                                  SF,
                                  typename XT::Common::multiplication_promotion<SF, RF>::type,
                                  r_r,
                                  r_rC,
                                  RF,
                                  RGV>(source_space, range_space, matrix.as_imp());
} // ... make_matrix_operator(...)


template <class GV, size_t r, size_t rC, class F, class M>
ConstMatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC, F>
make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space, const XT::LA::MatrixInterface<M>& matrix)
{
  return ConstMatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC, F>(
      space, space, matrix.as_imp());
}


/**
 * \brief Base class for linear operators which are assembled into a matrix.
 *
 * Similar to the GlobalAssembler, we derive from the XT::Grid::Walker and povide custom append() methods to allow to
 * add local element and intersection operators. In contrast to the GlobalAssembler we already hold the target matrix
 * (or create one of appropriate size and given sparsity pattern), into which we want to assemble. The operator is
 * assembled by walking over the given assembly_gid_view (which defaults to the one fom the given space). This allows to
 * assemble an operator only on a smaller grid view than the one given from the space (similar functionality could be
 * achieved by appending this operator to another walker and by providing an appropriate filter).
 *
 * \note One could achieve similar functionality by deriving from GlobalAssembler directly, which would slightly
 *       simplify the implementation of the append methods. However, we do not want to expose the other append methods
 *       of GlobalAssembler here (it should not be possible to append a local functional to an operator), but want to
 *       expose the ones from the XT::Grid::Walker (it should be possible to append other element or intersection
 *       operators or two-forms).
 *
 * \note See ConstMatrixBasedOperator and OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 * \sa ConstMatrixBasedOperator
 * \sa XT::Grid::Walker
 * \sa GlobalAssembler
 */
template <class M,
          class AssemblyGridView,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          class SGV = AssemblyGridView,
          class F = double,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RF = double,
          class RGV = SGV,
          class SV = typename XT::LA::Container<typename M::ScalarType, M::vector_type>::VectorType,
          class RV = SV>
class MatrixBasedOperator : XT::Common::StorageProvider<M>,
                            public ConstMatrixBasedOperator<M, SGV, s_r, s_rC, SF, F, r_r, r_rC, RF, RGV, SV, RV>,
                            public XT::Grid::Walker<AssemblyGridView>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");

  using ThisType = MatrixBasedOperator<M, AssemblyGridView, s_r, s_rC, SF, SGV, F, r_r, r_rC, RF, RGV, SV, RV>;
  using MatrixStorage = XT::Common::StorageProvider<M>;
  using OperatorBaseType = ConstMatrixBasedOperator<M, SGV, s_r, s_rC, SF, F, r_r, r_rC, RF, RGV, SV, RV>;
  using WalkerBaseType = XT::Grid::Walker<AssemblyGridView>;

public:
  using AssemblyGridViewType = AssemblyGridView;

  using typename OperatorBaseType::MatrixType;
  using typename OperatorBaseType::SourceSpaceType;
  using typename OperatorBaseType::RangeSpaceType;

  using typename WalkerBaseType::ElementType;
  using ElementFilterType = XT::Grid::ElementFilter<AssemblyGridViewType>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<AssemblyGridViewType>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<AssemblyGridViewType>;
  using ApplyOnAllIntersections = XT::Grid::ApplyOn::AllIntersections<AssemblyGridViewType>;

  using typename WalkerBaseType::E;
  using typename WalkerBaseType::I;

  /**
   * Ctor which accept an existing matrix into which to assemble.
   */
  MatrixBasedOperator(AssemblyGridViewType assembly_grid_view,
                      const SourceSpaceType& source_spc,
                      const RangeSpaceType& range_spc,
                      MatrixType& mat)
    : MatrixStorage(mat)
    , OperatorBaseType(source_spc, range_spc, MatrixStorage::access())
    , WalkerBaseType(assembly_grid_view)
    , assembled_(false)
  {
    // to detect assembly
    this->append([&](const auto&) { assembled_ = true; });
  }

  /**
   * Ctor which creates an appropriate matrix into which to assemble from a given sparsity pattern.
   */

  MatrixBasedOperator(AssemblyGridViewType assembly_grid_view,
                      const SourceSpaceType& source_spc,
                      const RangeSpaceType& range_spc,
                      const XT::LA::SparsityPatternDefault& pattern)
    : MatrixStorage(new MatrixType(range_spc.mapper().size(), source_spc.mapper().size(), pattern))
    , OperatorBaseType(source_spc, range_spc, MatrixStorage::access())
    , WalkerBaseType(assembly_grid_view)
    , assembled_(false)
  {
    // to detect assembly
    this->append([&](const auto&) { assembled_ = true; });
  }

  using OperatorBaseType::matrix;

  MatrixType& matrix()
  {
    return MatrixStorage::access();
  }

  using WalkerBaseType::append;

  ThisType& append(const LocalElementBilinearFormInterface<E, r_r, r_rC, RF, F, s_r, s_rC, SF>& local_bilinear_form,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    using LocalAssemblerType =
        LocalElementBilinearFormAssembler<MatrixType, AssemblyGridViewType, r_r, r_rC, RF, RGV, SGV, s_r, s_rC, SF>;
    this->append(new LocalAssemblerType(
                     this->range_space(), this->source_space(), local_bilinear_form, MatrixStorage::access(), param),
                 filter);
    return *this;
  }

  ThisType&
  append(const LocalIntersectionBilinearFormInterface<I, r_r, r_rC, RF, F, s_r, s_rC, SF>& local_bilinear_form,
         const XT::Common::Parameter& param = {},
         const IntersectionFilterType& filter = ApplyOnAllIntersections())
  {
    using LocalAssemblerType = LocalIntersectionBilinearFormAssembler<MatrixType,
                                                                      AssemblyGridViewType,
                                                                      r_r,
                                                                      r_rC,
                                                                      RF,
                                                                      RGV,
                                                                      SGV,
                                                                      s_r,
                                                                      s_rC,
                                                                      SF>;
    this->append(new LocalAssemblerType(
                     this->range_space(), this->source_space(), local_bilinear_form, MatrixStorage::access(), param),
                 filter);
    return *this;
  } // ... append(...)

  // similar append for LocalIntersectionFunctionalInterface ...

  void assemble(const bool use_tbb = false) override final
  {
    if (assembled_)
      return;
    // This clears all appended operators, which is ok, since we are done after assembling once!
    this->walk(use_tbb);
    assembled_ = true;
  }

private:
  bool assembled_;
}; // class MatrixBasedOperator


/// \name Variants of make_matrix_operator for a given matrix.
/// \{

template <class AGV,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class SF,
          class RGV,
          size_t r_r,
          size_t r_rC,
          class RF,
          class M>
MatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type,
                    GridView<AGV>,
                    s_r,
                    s_rC,
                    SF,
                    SGV,
                    typename XT::Common::multiplication_promotion<SF, RF>::type,
                    r_r,
                    r_rC,
                    RF,
                    RGV>
make_matrix_operator(GridView<AGV> assembly_grid_view,
                     const SpaceInterface<SGV, s_r, s_rC, SF>& source_space,
                     const SpaceInterface<RGV, r_r, r_rC, RF>& range_space,
                     XT::LA::MatrixInterface<M>& matrix)
{
  return MatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type,
                             GridView<AGV>,
                             s_r,
                             s_rC,
                             SF,
                             SGV,
                             typename XT::Common::multiplication_promotion<SF, RF>::type,
                             r_r,
                             r_rC,
                             RF,
                             RGV>(assembly_grid_view, source_space, range_space, matrix.as_imp());
} // ... make_matrix_operator(...)

template <class AGV, class GV, size_t r, size_t rC, class F, class M>
MatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type, GridView<AGV>, r, rC, F, GV>
make_matrix_operator(GridView<AGV> assembly_grid_view,
                     const SpaceInterface<GV, r, rC, F>& space,
                     XT::LA::MatrixInterface<M>& matrix)
{
  return MatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type, GridView<AGV>, r, rC, F, GV>(
      assembly_grid_view, space, space, matrix.as_imp());
}

template <class GV, size_t r, size_t rC, class F, class M>
MatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC, F>
make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space, XT::LA::MatrixInterface<M>& matrix)
{
  return MatrixBasedOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC, F>(
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
template <class MatrixType,
          class AGV,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class SF,
          class RGV,
          size_t r_r,
          size_t r_rC,
          class RF>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
                        MatrixBasedOperator<MatrixType,
                                            GridView<AGV>,
                                            s_r,
                                            s_rC,
                                            SF,
                                            SGV,
                                            typename XT::Common::multiplication_promotion<SF, RF>::type,
                                            r_r,
                                            r_rC,
                                            RF,
                                            RGV>>::type
make_matrix_operator(GridView<AGV> assembly_grid_view,
                     const SpaceInterface<SGV, s_r, s_rC, SF>& source_space,
                     const SpaceInterface<RGV, r_r, r_rC, RF>& range_space,
                     const XT::LA::SparsityPatternDefault& pattern)
{
  return MatrixBasedOperator<MatrixType,
                             GridView<AGV>,
                             s_r,
                             s_rC,
                             SF,
                             SGV,
                             typename XT::Common::multiplication_promotion<SF, RF>::type,
                             r_r,
                             r_rC,
                             RF,
                             RGV>(assembly_grid_view, source_space, range_space, pattern);
} // ... make_matrix_operator(...)

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(assembly_grid_view, space, pattern);
\endcode
 */
template <class MatrixType, class AGV, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
                        MatrixBasedOperator<MatrixType, GridView<AGV>, r, rC, F, GV>>::type
make_matrix_operator(GridView<AGV> assembly_grid_view,
                     const SpaceInterface<GV, r, rC, F>& space,
                     const XT::LA::SparsityPatternDefault& pattern)
{
  return MatrixBasedOperator<MatrixType, GridView<AGV>, r, rC, F, GV>(assembly_grid_view, space, space, pattern);
}

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(space, pattern);
\endcode
 */
template <class MatrixType, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value, MatrixBasedOperator<MatrixType, GV, r, rC, F>>::type
make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space, const XT::LA::SparsityPatternDefault& pattern)
{
  return MatrixBasedOperator<MatrixType, GV, r, rC, F>(space.grid_view(), space, space, pattern);
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
template <class MatrixType,
          class AGV,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class SF,
          class RGV,
          size_t r_r,
          size_t r_rC,
          class RF>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
                        MatrixBasedOperator<MatrixType,
                                            GridView<AGV>,
                                            s_r,
                                            s_rC,
                                            SF,
                                            SGV,
                                            typename XT::Common::multiplication_promotion<SF, RF>::type,
                                            r_r,
                                            r_rC,
                                            RF,
                                            RGV>>::type
make_matrix_operator(GridView<AGV> assembly_grid_view,
                     const SpaceInterface<SGV, s_r, s_rC, SF>& source_space,
                     const SpaceInterface<RGV, r_r, r_rC, RF>& range_space,
                     const Stencil stencil = Stencil::element_and_intersection)
{
  return MatrixBasedOperator<MatrixType,
                             GridView<AGV>,
                             s_r,
                             s_rC,
                             SF,
                             SGV,
                             typename XT::Common::multiplication_promotion<SF, RF>::type,
                             r_r,
                             r_rC,
                             RF,
                             RGV>(assembly_grid_view,
                                  source_space,
                                  range_space,
                                  make_sparsity_pattern(range_space, source_space, assembly_grid_view, stencil));
} // ... make_matrix_operator(...)

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(assembly_grid_view, space, stencil);
\endcode
 */
template <class MatrixType, class AGV, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
                        MatrixBasedOperator<MatrixType, GridView<AGV>, r, rC, F, GV>>::type
make_matrix_operator(GridView<AGV> assembly_grid_view,
                     const SpaceInterface<GV, r, rC, F>& space,
                     const Stencil stencil = Stencil::element_and_intersection)
{
  return MatrixBasedOperator<MatrixType, GridView<AGV>, r, rC, F, GV>(
      assembly_grid_view, space, space, make_sparsity_pattern(space, assembly_grid_view, stencil));
}

/**
 * \note Use as in
\code
auto op = make_matrix_operator<MatrixType>(space, stencil);
\endcode
 */
template <class MatrixType, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value, MatrixBasedOperator<MatrixType, GV, r, rC, F>>::type
make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space,
                     const Stencil stencil = Stencil::element_and_intersection)
{
  return MatrixBasedOperator<MatrixType, GV, r, rC, F>(
      space.grid_view(), space, space, make_sparsity_pattern(space, stencil));
}

/// \}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_MATRIX_BASED_OPERATOR_HH
