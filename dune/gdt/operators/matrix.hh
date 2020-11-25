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

#ifndef DUNE_GDT_OPERATORS_MATRIX_HH
#define DUNE_GDT_OPERATORS_MATRIX_HH

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
#include <dune/gdt/print.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


// forward, include is below
template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class SGV, class RGV>
class BilinearForm;


/**
 * \brief Base class for linear operators which are given by an assembled matrix.
 *
 * \note See OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 * \sa MatrixOperator
 */
template <class AGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class M = XT::LA::IstlRowMajorSparseMatrix<F>,
          class SGV = AGV,
          class RGV = AGV>
class ConstMatrixOperator : public OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>
{
public:
  using ThisType = ConstMatrixOperator;
  using BaseType = OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;

  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  ConstMatrixOperator(const AssemblyGridViewType& assembly_grid_vw,
                      const SourceSpaceType& source_spc,
                      const RangeSpaceType& range_spc,
                      const MatrixType& mat,
                      const std::string& logging_prefix = "",
                      const std::array<bool, 3>& logging_state = {{false, false, true}})
    : BaseType({}, logging_prefix.empty() ? "ConstMatrixOperator" : logging_prefix, logging_state)
    , assembly_grid_view_(assembly_grid_vw)
    , source_space_(source_spc)
    , range_space_(range_spc)
    , matrix_(mat)
    , linear_solver_(matrix_, source_space_.dof_communicator())
  {
    LOG_(debug) << "ConstMatrixOperator(source_space=" << &source_spc << ", range_space=" << &range_spc
                << ", matrix=" << &mat << ")" << std::endl;
    DUNE_THROW_IF(matrix_.rows() != range_space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "matrix_.rows() = " << matrix_.rows()
                                      << "\n   range_space_.mapper().size() = " << range_space_.mapper().size());
    DUNE_THROW_IF(matrix_.cols() != source_space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "matrix_.cols() = " << matrix_.cols()
                                      << "\n   source_space_.mapper().size() = " << source_space_.mapper().size());
  } // ConstMatrixOperator(...)

  ConstMatrixOperator(const ThisType& other)
    : BaseType(other)
    , assembly_grid_view_(other.assembly_grid_view_)
    , source_space_(other.source_space_)
    , range_space_(other.range_space_)
    , matrix_(other.matrix_)
    , linear_solver_(matrix_, source_space_.dof_communicator())
  {}

  ConstMatrixOperator(ThisType&& source)
    : BaseType(source)
    , assembly_grid_view_(source.assembly_grid_view_)
    , source_space_(source.source_space_)
    , range_space_(source.range_space_)
    , matrix_(source.matrix_)
    , linear_solver_(matrix_, source_space_.dof_communicator())
  {}

  // pull in methods from BilinearFormInterface, OperatorInterface, OperatorInterface
  using BaseType::apply;
  using BaseType::apply_inverse;
  using BaseType::jacobian;

  /// \name Required by OperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override
  {
    return range_space_;
  }

  bool linear() const override final
  {
    return true;
  }

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override
  {
    return source_space_;
  }

  const AssemblyGridViewType& assembly_grid_view() const override
  {
    return assembly_grid_view_;
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param = {}) const override
  {
    LOG_(debug) << "apply(source.sup_norm()=" << source.sup_norm() << ", range.sup_norm=" << range.sup_norm()
                << ", param=" << print(param, {{"oneline", "true"}}) << ")" << std::endl;
    this->assert_matching_source(source);
    this->assert_matching_range(range);
    try {
      LOG_(info) << "computing matrix*source ..." << std::endl;
      matrix_.mv(source, range);
    } catch (const XT::Common::Exceptions::shapes_do_not_match& ee) {
      DUNE_THROW(Exceptions::operator_error,
                 "when applying matrix to source and range dofs!\n\nThis was the original error: " << ee.what());
    }
  } // ... apply(...)

protected:
  std::vector<XT::Common::Configuration> all_jacobian_options() const override final
  {
    return {{{"type", "matrix"}}};
  }

public:
  void jacobian(const VectorType& source,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& param = {}) const override
  {
    LOG_(debug) << "jacobian(source.sup_norm()=" << source.sup_norm()
                << ", jacobian_op.matrix().sup_norm()=" << jacobian_op.matrix().sup_norm()
                << ", opts=" << print(opts, {{"oneline", "true"}}) << ", param=" << param << ")" << std::endl;
    this->assert_jacobian_opts(opts); // ensures that type matrix is requested
    LOG_(debug) << "   adding matrix_ * jacobian_op.scaling (matrix_.sup_norm() = " << matrix_.sup_norm()
                << ", jacobian_op.scaling = " << jacobian_op.scaling << ")" << std::endl;
    jacobian_op.matrix().axpy(jacobian_op.scaling, matrix_);
  } // ... jacobian(...)

protected:
  std::vector<XT::Common::Configuration> all_invert_options() const override final
  {
    std::vector<XT::Common::Configuration> ret;
    for (auto&& type : linear_solver_.types()) {
      auto opts = linear_solver_.options(type);
      opts["type"] = type;
      ret.emplace_back(std::move(opts));
    }
    for (auto&& opts : BaseType::all_invert_options())
      ret.emplace_back(opts);
    return ret;
  } // ... all_invert_options(...)

public:
  void apply_inverse(const VectorType& range,
                     VectorType& source,
                     const XT::Common::Configuration& opts,
                     const XT::Common::Parameter& param = {}) const override
  {
    LOG_(debug) << "apply_inverse(range.sup_norm()=" << range.sup_norm() << ", source.sup_norm()=" << source.sup_norm()
                << ", opts=" << print(opts, {{"oneline", "true"}}) << ", param=" << print(param, {{"oneline", "true"}})
                << ")" << std::endl;
    this->assert_apply_inverse_opts(opts);
    try {
      linear_solver_.apply(range, source, opts);
    } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
      DUNE_THROW(Exceptions::operator_error,
                 "when applying linear solver!\n\nThis was the original error: " << ee.what());
    } catch (const XT::Common::Exceptions::configuration_error&) {
      // we get this kind of error when type is not a linear solver type, it may thus be one from base
      BaseType::apply_inverse(range, source, opts, param);
    }
  } // ... apply_inverse(...)

  /// \}

  const MatrixType& matrix() const
  {
    return matrix_;
  }

protected:
  const AssemblyGridViewType& assembly_grid_view_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const MatrixType& matrix_;
  const XT::LA::Solver<MatrixType, typename SourceSpaceType::DofCommunicatorType> linear_solver_;
}; // class ConstMatrixOperator


template <class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC,
          class MatrixType>
auto make_matrix_operator(const AssemblyGridViewType& assembly_grid_view,
                          const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                          const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                          const MatrixType& matrix,
                          const std::string& logging_prefix = "")
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return ConstMatrixOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV>(
      assembly_grid_view, source_space, range_space, matrix, logging_prefix);
}


template <class GV, size_t r, size_t rC, class F, class MatrixType>
auto make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space,
                          const MatrixType& matrix,
                          const std::string& logging_prefix = "")
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return ConstMatrixOperator<GV, r, rC, r, rC, F, MatrixType, GV, GV>(
      space.grid_view(), space, space, matrix, logging_prefix);
}


/**
 * \brief Base class for linear operators which are assembled into a matrix.
 *
 * We derive from the XT::Grid::ElementAndIntersectionFunctor and povide custom append() methods to allow to add local
 * element and intersection operators. We also hold the target matrix into which we want to assemble. The operator is
 * assembled by walking over the given assembly_grid_view. If you want to assemble an operator only on a smaller grid
 * view, consider to append this operator to another walker or provide appropriate filters.
 *
 * \note See ConstMatrixOperator and OperatorInterface for a description of the template arguments.
 *
 * \todo Add logging to intersection and coupling intersection assemblers
 *
 * \todo Re-add derived logging for the local assemblers!
 * \todo Use std::pair instead of std::tuple!
 *
 * \note In compute_locally, the filters are evaluated w.r.t. assembly_grid_view and not the grid view
 *       of the walker this operator is appended to. This might not be what we want.
 *
 * \sa DisccreteOperatorInterface
 * \sa BilinearFormInterface
 * \sa ConstMatrixOperator
 * \sa XT::Grid::ElementAndIntersectionFunctor
 * \sa XT::Grid::Walker
 */
template <class AGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class M = XT::LA::IstlRowMajorSparseMatrix<F>,
          class SGV = AGV,
          class RGV = AGV>
class MatrixOperator
  : XT::Common::StorageProvider<M>
  , public ConstMatrixOperator<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>
  , public XT::Grid::ElementAndIntersectionFunctor<AGV>
{
  using MatrixStorage = XT::Common::StorageProvider<M>;

public:
  using ThisType = MatrixOperator;
  using BaseOperatorType = ConstMatrixOperator<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;
  using BaseFunctorType = XT::Grid::ElementAndIntersectionFunctor<AGV>;

  using BaseOperatorType::logger;
  using typename BaseOperatorType::AssemblyGridViewType;
  using typename BaseOperatorType::FieldType;
  using typename BaseOperatorType::MatrixType;
  using typename BaseOperatorType::RangeSpaceType;
  using typename BaseOperatorType::SourceSpaceType;
  using typename BaseOperatorType::V;
  using typename BaseOperatorType::VectorType;

  using typename BaseFunctorType::E;
  using typename BaseFunctorType::I;

  using ElementFilterType = XT::Grid::ElementFilter<AGV>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<AGV>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<AGV>;
  using ApplyOnAllIntersections = XT::Grid::ApplyOn::AllIntersections<AGV>;

protected:
  using LocalElementBilinearFormAssemblerType =
      LocalElementBilinearFormAssembler<M, AGV, r_r, r_rC, F, RGV, SGV, s_r, s_rC, F>;
  using LocalCouplingIntersectionBilinearFormAssemblerType =
      LocalCouplingIntersectionBilinearFormAssembler<M, AGV, r_r, r_rC, F, RGV, SGV, s_r, s_rC, F>;
  using LocalIntersectionBilinearFormAssemblerType =
      LocalIntersectionBilinearFormAssembler<M, AGV, r_r, r_rC, F, RGV, SGV, s_r, s_rC, F>;

  using LocalElementOperatorFiniteDifferenceJacobianAssemblerType =
      LocalElementOperatorFiniteDifferenceJacobianAssembler<M, SGV, s_r, s_rC, F, r_r, r_rC, RGV>;
  using LocalIntersectionOperatorFiniteDifferenceJacobianAssemblerType =
      LocalIntersectionOperatorFiniteDifferenceJacobianAssembler<M, SGV, s_r, s_rC, F, r_r, r_rC, RGV>;

public:
  using BilinearFormType = BilinearForm<AGV, s_r, s_rC, r_r, r_rC, F, SGV, RGV>;

  /// \brief Defines the scaling of the to-be-appended jacobian contributions
  /// \sa    ConstantMatrixOperator::jacobian
  /// \sa    OperatorInterface::apply_inverse
  FieldType scaling;

  MatrixOperator(const AssemblyGridViewType& assembly_grid_vw,
                 const SourceSpaceType& source_spc,
                 const RangeSpaceType& range_spc,
                 MatrixType& mat,
                 const std::string& logging_prefix = "",
                 const std::array<bool, 3>& logging_state = {{false, false, true}})
    : MatrixStorage(mat)
    , BaseOperatorType(assembly_grid_vw,
                       source_spc,
                       range_spc,
                       MatrixStorage::access(),
                       logging_prefix.empty() ? "MatrixOperator" : logging_prefix,
                       logging_state)
    , BaseFunctorType(logging_prefix.empty() ? "MatrixOperator" : logging_prefix, logging_state)
    , scaling(1.)
  {
    LOG_(debug) << "MatrixOperator(assembly_grid_view=" << &assembly_grid_vw << ", source_space=" << &source_spc
                << ", range_space=" << &range_spc << ", mat=" << &mat << ")" << std::endl;
  }

  MatrixOperator(const AssemblyGridViewType& assembly_grid_vw,
                 const SourceSpaceType& source_spc,
                 const RangeSpaceType& range_spc,
                 MatrixType*&& mat_ptr,
                 const std::string& logging_prefix = "",
                 const std::array<bool, 3>& logging_state = {{false, false, true}})
    : MatrixStorage(std::move(mat_ptr))
    , BaseOperatorType(assembly_grid_vw,
                       source_spc,
                       range_spc,
                       MatrixStorage::access(),
                       logging_prefix.empty() ? "MatrixOperator" : logging_prefix,
                       logging_state)
    , BaseFunctorType(logging_prefix.empty() ? "MatrixOperator" : logging_prefix, logging_state)
    , scaling(1.)
  {
    LOG_(debug) << "MatrixOperator(assembly_grid_view=" << &assembly_grid_vw << ", source_space=" << &source_spc
                << ", range_space=" << &range_spc << ", mat_ptr=" << &mat_ptr << ")" << std::endl;
  }

  MatrixOperator(const ThisType&) = delete;

  /// \note Performs something like a shallow copy (as required by ElementAndIntersectionFunctor), i.e. the copied
  ///       operator shares the matrix.
  /// \note Manual copy required to duplicate the local assemblers
  MatrixOperator(ThisType& other)
    : MatrixStorage(other)
    , BaseOperatorType(other)
    , BaseFunctorType(other)
    , scaling(other.scaling)
  {
    const auto copy_local_data = [](const auto& origin, auto& target) {
      for (auto& data : origin) {
        auto& local_assembler = *std::get<0>(data);
        const auto& filter = *std::get<1>(data);
        using LAT = std::decay_t<decltype(local_assembler)>;
        target.emplace_back(new LAT(local_assembler), filter.copy());
      }
    };
    copy_local_data(other.element_bilinear_form_data_, element_bilinear_form_data_);
    copy_local_data(other.element_fd_operator_data_, element_fd_operator_data_);
    copy_local_data(other.coupling_intersection_bilinear_form_data_, coupling_intersection_bilinear_form_data_);
    copy_local_data(other.intersection_bilinear_form_data_, intersection_bilinear_form_data_);
    copy_local_data(other.intersection_fd_operator_data_, intersection_fd_operator_data_);
  } // ... MatrixOperator(...)

  MatrixOperator(ThisType&& source) = default;

  // pull in methods from various base classes
  using BaseOperatorType::matrix;

  /// \name Required by ElementAndIntersectionFunctor.
  /// \{

  BaseFunctorType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const E& element) override final
  {
    const auto apply_local_assembler = [&](const auto& element_data) {
      for (auto& data : element_data) {
        auto& local_assembler = *std::get<0>(data);
        const auto& filter = *std::get<1>(data);
        if (filter.contains(this->assembly_grid_view_, element))
          local_assembler.apply_local(element);
      }
    };
    apply_local_assembler(element_bilinear_form_data_);
    apply_local_assembler(element_fd_operator_data_);
  } // ... apply_local(...)

  void apply_local(const I& intersection, const E& inside_element, const E& outside_element) override final
  {
    const auto apply_local_assembler = [&](auto& intersection_data) {
      for (auto& data : intersection_data) {
        auto& local_assembler = *std::get<0>(data);
        const auto& filter = *std::get<1>(data);
        if (filter.contains(this->assembly_grid_view_, intersection))
          local_assembler.apply_local(intersection, inside_element, outside_element);
      }
    };
    apply_local_assembler(coupling_intersection_bilinear_form_data_);
    apply_local_assembler(intersection_bilinear_form_data_);
    apply_local_assembler(intersection_fd_operator_data_);
  } // ... apply_local(...)

  void finalize() override final
  {
    // clear everything after assembling into the matrix
    element_bilinear_form_data_.clear();
    coupling_intersection_bilinear_form_data_.clear();
    intersection_bilinear_form_data_.clear();
    element_fd_operator_data_.clear();
    intersection_fd_operator_data_.clear();
  } // ... finalize(...)

  /// \}
  /// \name Required by BilinearFormInterface
  /// \{

  void assemble(const bool use_tbb = false) override final
  {
    XT::Grid::Walker<AGV> walker(this->assembly_grid_view_);
    walker.append(*this);
    walker.walk(use_tbb);
  }

  /// \}
  /// \name These methods extend the ones from ConstantMatrixOperator
  /// \{

  MatrixType& matrix()
  {
    return MatrixStorage::access();
  }

  /// \}
  /// \name These methods allow to append the local bilinear forms from a BilinearForm
  /// \{

  ThisType& append(const BilinearFormType& bilinear_form, const XT::Common::Parameter& param = {})
  {
    for (auto& data : bilinear_form.element_data()) {
      const auto& local_bilinear_form = *std::get<0>(data);
      const auto& filter = *std::get<1>(data);
      element_bilinear_form_data_.emplace_back(
          new LocalElementBilinearFormAssemblerType(this->range_space(),
                                                    this->source_space(),
                                                    local_bilinear_form,
                                                    MatrixStorage::access(),
                                                    param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
          filter.copy());
    }
    for (auto& data : bilinear_form.coupling_intersection_data()) {
      const auto& local_bilinear_form = *std::get<0>(data);
      const auto& filter = *std::get<1>(data);
      coupling_intersection_bilinear_form_data_.emplace_back(
          new LocalCouplingIntersectionBilinearFormAssemblerType(
              this->range_space(),
              this->source_space(),
              local_bilinear_form,
              MatrixStorage::access(),
              param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
          filter.copy());
    }
    for (auto& data : bilinear_form.intersection_data()) {
      const auto& local_bilinear_form = *std::get<0>(data);
      const auto& filter = *std::get<1>(data);
      intersection_bilinear_form_data_.emplace_back(
          new LocalIntersectionBilinearFormAssemblerType(
              this->range_space(),
              this->source_space(),
              local_bilinear_form,
              MatrixStorage::access(),
              param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
          filter.copy());
    }
    return *this;
  } // ... append(...)

  /// \name These methods allow to compute the jacobian of the appended local operators by finite differences.
  /// \{

  ThisType& append(const LocalElementOperatorInterface<V, SGV, s_r, s_rC, F, r_r, r_rC, F, RGV>& local_operator,
                   const VectorType& source,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    element_fd_operator_data_.emplace_back(new LocalElementOperatorFiniteDifferenceJacobianAssemblerType(
                                               this->source_space(),
                                               this->range_space(),
                                               MatrixStorage::access(),
                                               source,
                                               local_operator,
                                               param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
                                           filter.copy());
    return *this;
  } // ... append(...)

  ThisType& append(const LocalIntersectionOperatorInterface<I, V, SGV, s_r, s_rC, F, r_r, r_rC>& local_operator,
                   const VectorType& source,
                   const XT::Common::Parameter& param = {},
                   const IntersectionFilterType& filter = ApplyOnAllIntersections())
  {
    intersection_fd_operator_data_.emplace_back(new LocalIntersectionOperatorFiniteDifferenceJacobianAssemblerType(
                                                    this->source_space(),
                                                    this->range_space(),
                                                    MatrixStorage::access(),
                                                    source,
                                                    local_operator,
                                                    param + XT::Common::Parameter("matrixoperator.scaling", scaling)),
                                                filter.copy());
    return *this;
  } // ... append(...)

  /// \}

protected:
  std::list<std::tuple<std::unique_ptr<LocalElementBilinearFormAssemblerType>, std::unique_ptr<ElementFilterType>>>
      element_bilinear_form_data_;
  std::list<std::tuple<std::unique_ptr<LocalCouplingIntersectionBilinearFormAssemblerType>,
                       std::unique_ptr<IntersectionFilterType>>>
      coupling_intersection_bilinear_form_data_;
  std::list<
      std::tuple<std::unique_ptr<LocalIntersectionBilinearFormAssemblerType>, std::unique_ptr<IntersectionFilterType>>>
      intersection_bilinear_form_data_;
  std::list<std::tuple<std::unique_ptr<LocalElementOperatorFiniteDifferenceJacobianAssemblerType>,
                       std::unique_ptr<ElementFilterType>>>
      element_fd_operator_data_;
  std::list<std::tuple<std::unique_ptr<LocalIntersectionOperatorFiniteDifferenceJacobianAssemblerType>,
                       std::unique_ptr<IntersectionFilterType>>>
      intersection_fd_operator_data_;
}; // class MatrixOperator


/// \name Variants of make_matrix_operator for a given matrix.
/// \{


template <class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC,
          class MatrixType>
auto make_matrix_operator(const AssemblyGridViewType& assembly_grid_view,
                          const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                          const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                          MatrixType& matrix,
                          const std::string& logging_prefix = "")
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return MatrixOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV>(
      assembly_grid_view, source_space, range_space, matrix, logging_prefix);
}


template <class GV, size_t r, size_t rC, class F, class MatrixType>
auto make_matrix_operator(const SpaceInterface<GV, r, rC, F>& space,
                          MatrixType& matrix,
                          const std::string& logging_prefix = "")
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return MatrixOperator<GV, r, rC, r, rC, F, MatrixType, GV, GV>(
      space.grid_view(), space, space, matrix, logging_prefix);
}


/// \}
/// \name Variants of make_matrix_operator, where an appropriate matrix is created from a given pattern
/// \{


template <class AssemblyGridViewType, class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC>
auto make_matrix_operator(const AssemblyGridViewType& assembly_grid_view,
                          const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                          const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                          const XT::LA::SparsityPatternDefault& pattern,
                          const std::string& logging_prefix = "")
{
  using M = XT::LA::IstlRowMajorSparseMatrix<F>;
  return MatrixOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>(
      assembly_grid_view,
      source_space,
      range_space,
      new M(range_space.mapper().size(), source_space.mapper().size(), pattern),
      logging_prefix);
} // ... make_matrix_operator(...)


/**
\note Use as in
\code
auto op = make_matrix_operator<MatrixType>(assembly_grid_view, source_space, range_space, pattern);
\endcode
 */
template <class MatrixType,
          class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC>
auto make_matrix_operator(const AssemblyGridViewType& assembly_grid_view,
                          const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                          const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                          const XT::LA::SparsityPatternDefault& pattern,
                          const std::string& logging_prefix = "")
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return MatrixOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV>(
      assembly_grid_view,
      source_space,
      range_space,
      new MatrixType(range_space.mapper().size(), source_space.mapper().size(), pattern),
      logging_prefix);
} // ... make_matrix_operator(...)


template <class GV, size_t s_r, size_t s_rC, class F, size_t r_r, size_t r_rC>
auto make_matrix_operator(const SpaceInterface<GV, r_r, r_rC, F>& space,
                          const XT::LA::SparsityPatternDefault& pattern,
                          const std::string& logging_prefix = "")
{
  using M = XT::LA::IstlRowMajorSparseMatrix<F>;
  return MatrixOperator<GV, s_r, s_rC, r_r, r_rC, F, M, GV, GV>(
      space.grid_view(), space, space, new M(space.mapper().size(), space.mapper().size(), pattern), logging_prefix);
}


/**
 * \note Use as in
\code
 auto op = make_matrix_operator<MatrixType>(space, pattern);
\endcode
 */
template <class MatrixType, class GV, size_t s_r, size_t s_rC, class F, size_t r_r, size_t r_rC>
auto make_matrix_operator(const SpaceInterface<GV, r_r, r_rC, F>& space,
                          const XT::LA::SparsityPatternDefault& pattern,
                          const std::string& logging_prefix = "")
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return MatrixOperator<GV, s_r, s_rC, r_r, r_rC, F, MatrixType, GV, GV>(
      space.grid_view(),
      space,
      space,
      new MatrixType(space.mapper().size(), space.mapper().size(), pattern),
      logging_prefix);
} // ... make_matrix_operator(...)


/// \}
/// \name Variants of make_matrix_operator, where an appropriate matrix is created from given stencil
/// \{


template <class AssemblyGridViewType, class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC>
auto make_matrix_operator(const AssemblyGridViewType& assembly_grid_view,
                          const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                          const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                          const Stencil stencil = Stencil::element_and_intersection,
                          const std::string& logging_prefix = "")
{
  using M = XT::LA::IstlRowMajorSparseMatrix<F>;
  return MatrixOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>(
      assembly_grid_view,
      source_space,
      range_space,
      new M(range_space.mapper().size(),
            source_space.mapper().size(),
            make_sparsity_pattern(range_space, source_space, assembly_grid_view, stencil)),
      logging_prefix);
} // ... make_matrix_operator(...)


/**
 * \note Use as in
\code
 auto op = make_matrix_operator<MatrixType>(assembly_grid_view, source_space, range_space, stencil);
\endcode
 */
template <class MatrixType,
          class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC>
auto make_matrix_operator(const AssemblyGridViewType& assembly_grid_view,
                          const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                          const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                          const Stencil stencil = Stencil::element_and_intersection,
                          const std::string& logging_prefix = "")
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return MatrixOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV>(
      assembly_grid_view,
      source_space,
      range_space,
      new MatrixType(range_space.mapper().size(),
                     source_space.mapper().size(),
                     make_sparsity_pattern(range_space, source_space, assembly_grid_view, stencil)),
      logging_prefix);
} // ... make_matrix_operator(...)


template <class GV, size_t s_r, size_t s_rC, class F, size_t r_r, size_t r_rC>
auto make_matrix_operator(const SpaceInterface<GV, r_r, r_rC, F>& space,
                          const Stencil stencil = Stencil::element_and_intersection,
                          const std::string& logging_prefix = "")
{
  using M = XT::LA::IstlRowMajorSparseMatrix<F>;
  return MatrixOperator<GV, s_r, s_rC, r_r, r_rC, F, M, GV, GV>(
      space.grid_view(),
      space,
      space,
      new M(space.mapper().size(), space.mapper().size(), make_sparsity_pattern(space, stencil)),
      logging_prefix);
} // ... make_matrix_operator(...)


/**
 * \note Use as in
\code
 auto op = make_matrix_operator<MatrixType>(space, stencil);
\endcode
 */
template <class MatrixType, class GV, size_t s_r, size_t s_rC, class F, size_t r_r, size_t r_rC>
auto make_matrix_operator(const SpaceInterface<GV, r_r, r_rC, F>& space,
                          const Stencil stencil = Stencil::element_and_intersection,
                          const std::string& logging_prefix = "")
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return MatrixOperator<GV, s_r, s_rC, r_r, r_rC, F, MatrixType, GV, GV>(
      space.grid_view(),
      space,
      space,
      new MatrixType(space.mapper().size(), space.mapper().size(), make_sparsity_pattern(space, stencil)),
      logging_prefix);
} // ... make_matrix_operator(...)


/// \}


} // namespace GDT
} // namespace Dune

#include "bilinear-form.hh"

#endif // DUNE_GDT_OPERATORS_MATRIX_HH
