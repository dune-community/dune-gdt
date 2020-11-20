// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_OPERATORS_BILINEAR_FORM_HH
#define DUNE_GDT_OPERATORS_BILINEAR_FORM_HH

#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/functions/grid-function.hh>

#include <dune/gdt/local/bilinear-forms/interfaces.hh>
#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {


// forwards
template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class SGV, class RGV>
class BilinearFormAssembler;

template <class AGV,
          size_t s_r,
          size_t s_rC,
          size_t r_r,
          size_t r_rC,
          class F,
          class M,
          class SGV,
          class RGV>
class MatrixOperator; // include is below


/**
 * \note We could also allow operator+=() to optionally take a Parameter, to fix the parameter value for this appended
 *       local bilinear form. We would need to extend the list members, not add this parameter dependency and check
 *       in the BilinearFormAssembler for each bilinear form is a parater was provided on operator+=() or if the
 *       parameter provided in apply2() is to be used. We do not yet do so, since the current approach is simpler.
 * \todo Use std::pair instead of std::tuple!
 */
template <class AssemblyGridView,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class SGV = AssemblyGridView,
          class RGV = AssemblyGridView>
class BilinearForm : public BilinearFormInterface<SGV, s_r, s_rC, r_r, r_rC, F, RGV>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");

  using ThisType = BilinearForm;
  using BaseType = BilinearFormInterface<SGV, s_r, s_rC, r_r, r_rC, F, RGV>;

public:
  using AssemblyGridViewType = AssemblyGridView;
  using AGV = AssemblyGridView;
  using E = XT::Grid::extract_entity_t<AGV>;
  using I = XT::Grid::extract_intersection_t<AGV>;
  using ElementFilterType = XT::Grid::ElementFilter<AGV>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<AGV>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<AGV>;
  using ApplyOnAllIntersections = XT::Grid::ApplyOn::AllIntersections<AGV>;

  using typename BaseType::FieldType;
  using typename BaseType::RangeFunctionType;
  using typename BaseType::SourceFunctionType;

  using LocalElementBilinearFormType = LocalElementBilinearFormInterface<E, r_r, r_rC, F, F, s_r, s_rC, F>;
  using LocalCouplingIntersectionBilinearFormType =
      LocalCouplingIntersectionBilinearFormInterface<I, r_r, r_rC, F, F, s_r, s_rC, F>;
  using LocalIntersectionBilinearFormType = LocalIntersectionBilinearFormInterface<I, r_r, r_rC, F, F, s_r, s_rC, F>;

  BilinearForm(const AssemblyGridViewType& assembly_grid_view,
               const std::string& logging_prefix = "",
               const std::array<bool, 3>& logging_state = {false, false, true})
    : BaseType({}, logging_prefix.empty() ? "BilinearForm" : logging_prefix, logging_state)
    , assembly_grid_view_(assembly_grid_view)
  {
    LOG_(debug) << "BilinearForm(assembly_grid_view=" << &assembly_grid_view << ")" << std::endl;
  }

  BilinearForm(const ThisType& other)
    : BaseType(other)
    , assembly_grid_view_(other.assembly_grid_view_)
  {
    const auto copy_local_data = [](const auto& origin, auto& target) {
      for (const auto& data : origin) {
        const auto& local_bilinear_form = *std::get<0>(data);
        const auto& filter = *std::get<1>(data);
        target.emplace_back(local_bilinear_form.copy(), filter.copy());
      }
    };
    copy_local_data(other.element_data_, element_data_);
    copy_local_data(other.coupling_intersection_data_, coupling_intersection_data_);
    copy_local_data(other.intersection_data_, intersection_data_);
  } // BilinearForm(...)

  BilinearForm(ThisType&&) = default;

  /// \brief allows to fix the arguments to apply2(), the resulting assembler can be appended to a GridWalker
  auto with(RangeFunctionType range_function,
            SourceFunctionType source_function,
            const XT::Common::Parameter& param = {}) const
  {
    return BilinearFormAssembler<AGV, s_r, s_rC, r_r, r_rC, F, SGV, RGV>(
        *this, range_function, source_function, param, this->logger.prefix + "_assembler", this->logger.state);
  }

  /// \brief creates the associated MatrixOperator (which still has to be assembled)
  template <class MatrixType>
  auto with(const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
            const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
            MatrixType& matrix,
            const XT::Common::Parameter& param = {}) const
  {
    static_assert(XT::LA::is_matrix<MatrixType>::value, "");
    MatrixOperator<AGV, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV> op(
        assembly_grid_view_, source_space, range_space, matrix, this->logger.prefix + "_mat", this->logger.state);
    op.append(*this, param);
    return op;
  } // ... with(...)

  /// \brief creates the associated MatrixOperator (which still has to be assembled)
  template <class MatrixType = XT::LA::IstlRowMajorSparseMatrix<F>>
  auto with(const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
            const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
            const XT::Common::Parameter& param = {}) const
  {
    static_assert(XT::LA::is_matrix<MatrixType>::value, "");
    MatrixOperator<AGV, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV> op(
        assembly_grid_view_,
        source_space,
        range_space,
        new MatrixType(range_space.mapper().size(),
                       source_space.mapper().size(),
                       make_element_and_intersection_sparsity_pattern(range_space, source_space, assembly_grid_view_)),
        this->logger.prefix + "_mat",
        this->logger.state);
    op.append(*this, param);
    return op;
  } // ... with(...)

  /// \name Required by BilinearFormInterface.
  /// \{

  FieldType apply2(RangeFunctionType range_function,
                   SourceFunctionType source_function,
                   const XT::Common::Parameter& param = {}) const override final
  {
    auto assembler = this->with(range_function, source_function, param);
    XT::Grid::Walker<AGV> walker(this->assembly_grid_view_);
    walker.append(assembler);
    walker.walk(/*use_tbb=*/true);
    return assembler.result();
  }

  /// \}

  const AssemblyGridViewType& assembly_grid_view() const
  {
    return assembly_grid_view_;
  }

  /// \name These methods allow to append local element bilinear forms
  ///\{

  ThisType& operator+=(const LocalElementBilinearFormType& local_bilinear_form)
  {
    LOG_(info) << "+=(local_element_bilinear_form=" << &local_bilinear_form << ")" << std::endl;
    this->extend_parameter_type(local_bilinear_form.parameter_type());
    element_data_.emplace_back(local_bilinear_form.copy(), new ApplyOnAllElements());
    return *this;
  }

  ThisType& operator+=(
      const std::tuple<const LocalElementBilinearFormType&, const ElementFilterType&>& local_bilinear_form__filter)
  {
    const auto& local_bilinear_form = std::get<0>(local_bilinear_form__filter);
    const auto& filter = std::get<1>(local_bilinear_form__filter);
    LOG_(info) << "+=(local_element_bilinear_form=" << &local_bilinear_form << ", filter=" << &filter << ")"
               << std::endl;
    this->extend_parameter_type(local_bilinear_form.parameter_type());
    element_data_.emplace_back(local_bilinear_form.copy(), filter.copy());
    return *this;
  } // ... operator+=(...)

  ///\}
  /// \name These methods allow to append local coupling intersection bilinear forms
  ///\{

  ThisType& operator+=(const LocalCouplingIntersectionBilinearFormType& local_bilinear_form)
  {
    LOG_(info) << "+=(local_coupling_intersection_bilinear_form=" << &local_bilinear_form << ")" << std::endl;
    this->extend_parameter_type(local_bilinear_form.parameter_type());
    coupling_intersection_data_.emplace_back(local_bilinear_form.copy(), new ApplyOnAllIntersections());
    return *this;
  }

  ThisType& operator+=(const std::tuple<const LocalCouplingIntersectionBilinearFormType&,
                                        const IntersectionFilterType&>& local_bilinear_form__filter)
  {
    const auto& local_bilinear_form = std::get<0>(local_bilinear_form__filter);
    const auto& filter = std::get<1>(local_bilinear_form__filter);
    LOG_(info) << "+=(local_coupling_intersection_bilinear_form=" << &local_bilinear_form << ", filter=" << &filter
               << ")" << std::endl;
    this->extend_parameter_type(local_bilinear_form.parameter_type());
    coupling_intersection_data_.emplace_back(local_bilinear_form.copy(), filter.copy());
    return *this;
  } // ... operator+=(...)

  ///\}
  /// \name These methods allow to append local intersection bilinear forms
  ///\{

  ThisType& operator+=(const LocalIntersectionBilinearFormType& local_bilinear_form)
  {
    LOG_(info) << "+=(local_intersection_bilinear_form=" << &local_bilinear_form << ")" << std::endl;
    this->extend_parameter_type(local_bilinear_form.parameter_type());
    intersection_data_.emplace_back(local_bilinear_form.copy(), new ApplyOnAllIntersections());
    return *this;
  }

  ThisType& operator+=(const std::tuple<const LocalIntersectionBilinearFormType&, const IntersectionFilterType&>&
                           local_bilinear_form__filter)
  {
    const auto& local_bilinear_form = std::get<0>(local_bilinear_form__filter);
    const auto& filter = std::get<1>(local_bilinear_form__filter);
    LOG_(info) << "+=(local_intersection_bilinear_form=" << &local_bilinear_form << ", filter=" << &filter << ")"
               << std::endl;
    this->extend_parameter_type(local_bilinear_form.parameter_type());
    intersection_data_.emplace_back(local_bilinear_form.copy(), filter.copy());
    return *this;
  }

  ///\}

  const std::list<std::tuple<std::unique_ptr<LocalElementBilinearFormType>, std::unique_ptr<ElementFilterType>>>&
  element_data() const
  {
    return element_data_;
  }

  const std::list<
      std::tuple<std::unique_ptr<LocalCouplingIntersectionBilinearFormType>, std::unique_ptr<IntersectionFilterType>>>&
  coupling_intersection_data() const
  {
    return coupling_intersection_data_;
  }

  const std::list<
      std::tuple<std::unique_ptr<LocalIntersectionBilinearFormType>, std::unique_ptr<IntersectionFilterType>>>&
  intersection_data() const
  {
    return intersection_data_;
  }

protected:
  const AssemblyGridViewType& assembly_grid_view_;
  std::list<std::tuple<std::unique_ptr<LocalElementBilinearFormType>, std::unique_ptr<ElementFilterType>>>
      element_data_;
  std::list<
      std::tuple<std::unique_ptr<LocalCouplingIntersectionBilinearFormType>, std::unique_ptr<IntersectionFilterType>>>
      coupling_intersection_data_;
  std::list<std::tuple<std::unique_ptr<LocalIntersectionBilinearFormType>, std::unique_ptr<IntersectionFilterType>>>
      intersection_data_;
}; // class BilinearForm


/**
 * \note In compute_locally, the filters are evaluated w.r.t. bilinear_form_.assembly_grid_view_ and not the grid view
 *       of the walker this assembler is appended to. This might not be what we want.
 */
template <class AGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class SGV = AGV,
          class RGV = AGV>
class BilinearFormAssembler
  : public XT::Grid::ElementAndIntersectionFunctor<AGV>
  , XT::Common::ThreadResultPropagator<BilinearFormAssembler<AGV, s_r, s_rC, r_r, r_rC, F, SGV, RGV>, F>
{
  using ThisType = BilinearFormAssembler;
  using BaseType = XT::Grid::ElementAndIntersectionFunctor<AGV>;
  using Propagator =
      XT::Common::ThreadResultPropagator<BilinearFormAssembler<AGV, s_r, s_rC, r_r, r_rC, F, SGV, RGV>, F>;
  friend Propagator;

public:
  using BilinearFormType = BilinearForm<AGV, s_r, s_rC, r_r, r_rC, F, SGV, RGV>;
  using E = typename BilinearFormType::E;
  using I = typename BilinearFormType::I;
  using FieldType = typename BilinearFormType::FieldType;
  using SourceFunctionType = typename BilinearFormType::SourceFunctionType;
  using RangeFunctionType = typename BilinearFormType::RangeFunctionType;

  BilinearFormAssembler(BilinearFormType bilinear_form,
                        RangeFunctionType rng,
                        SourceFunctionType src,
                        const XT::Common::Parameter& param = {},
                        const std::string& logging_prefix = "",
                        const std::array<bool, 3>& logging_state = {false, false, true})
    : BaseType(logging_prefix.empty() ? "BilinearFormAssembler" : logging_prefix, logging_state)
    , Propagator(this)
    , bilinear_form_(bilinear_form)
    , range_(rng.copy_as_grid_function())
    , source_(src.copy_as_grid_function())
    , param_(param)
    , local_range_in_(range_->local_function())
    , local_range_out_(range_->local_function())
    , local_source_in_(source_->local_function())
    , local_source_out_(source_->local_function())
    , bilinear_form_value_in_in_(1, 1, 0.)
    , bilinear_form_value_in_out_(1, 1, 0.)
    , bilinear_form_value_out_in_(1, 1, 0.)
    , bilinear_form_value_out_out_(1, 1, 0.)
    , result_(0)
  {}

  BilinearFormAssembler(const ThisType& other)
    : BaseType(other)
    , Propagator(this)
    , bilinear_form_(other.bilinear_form_)
    , range_(other.range_->copy_as_grid_function())
    , source_(other.source_->copy_as_grid_function())
    , param_(other.param_)
    , local_range_in_(range_->local_function())
    , local_range_out_(range_->local_function())
    , local_source_in_(source_->local_function())
    , local_source_out_(source_->local_function())
    , bilinear_form_value_in_in_(other.bilinear_form_value_in_in_)
    , bilinear_form_value_in_out_(other.bilinear_form_value_in_out_)
    , bilinear_form_value_out_in_(other.bilinear_form_value_out_in_)
    , bilinear_form_value_out_out_(other.bilinear_form_value_out_out_)
    , result_(other.result_)
  {}

  BilinearFormAssembler(ThisType&& source) = default;

  /// \name Required by ElementAndIntersectionFunctor.
  /// \{

  BaseType* copy() override final
  {
    return Propagator::copy_imp();
  }

  void prepare() override final
  {
    result_ = 0;
  }

  void apply_local(const E& element) override final
  {
    result_ += compute_locally(element);
  }

  void apply_local(const I& intersection, const E& inside_element, const E& outside_element) override final
  {
    result_ += compute_locally(intersection, inside_element, outside_element);
  }

  void finalize() override final
  {
    Propagator::finalize_imp();
  }

  /// \}
  /// \name Required by ThreadResultPropagator.
  /// \{

protected:
  void set_result(FieldType res)
  {
    result_ = res;
  }

public:
  /// \}
  /// \name Can be used to compute the application of the bilinear form to source and range locally.
  /// \{

  /// \brief Variant of compute_locally to apply all local element bilinear forms
  FieldType compute_locally(const E& element)
  {
    FieldType ret = 0;
    local_source_in_->bind(element);
    local_range_in_->bind(element);
    for (auto& data : bilinear_form_.element_data()) {
      const auto& local_bilinear_form = *std::get<0>(data);
      const auto& filter = *std::get<1>(data);
      if (filter.contains(bilinear_form_.assembly_grid_view(), element)) {
        assert(local_source_in_->size(param_) == 1);
        assert(local_range_in_->size(param_) == 1);
        local_bilinear_form.apply2(*local_range_in_, *local_source_in_, bilinear_form_value_in_in_, param_);
        assert(bilinear_form_value_in_in_.rows() >= 1);
        assert(bilinear_form_value_in_in_.cols() >= 1);
        ret += bilinear_form_value_in_in_[0][0];
      }
    }
    return ret;
  } // ... compute_locally(...)

  /// \brief Variant of compute_locally to apply all local intersection and coupling intersection bilinear forms
  ///
  /// \attention Make sure to call this method only on suitable intersections, neither intersection.neighbor() not
  ///            the validity of outside_element are checked!
  FieldType compute_locally(const I& intersection, const E& inside_element, const E& outside_element)
  {
    FieldType ret = 0;
    local_source_in_->bind(inside_element);
    local_range_in_->bind(inside_element);
    local_source_out_->bind(outside_element);
    local_range_out_->bind(outside_element);
    for (auto& data : bilinear_form_.coupling_intersection_data()) {
      const auto& local_bilinear_form = *std::get<0>(data);
      const auto& filter = *std::get<1>(data);
      if (filter.contains(bilinear_form_.assembly_grid_view(), intersection)) {
        assert(local_source_in_->size(param_) == 1);
        assert(local_source_out_->size(param_) == 1);
        assert(local_range_in_->size(param_) == 1);
        assert(local_range_out_->size(param_) == 1);
        local_bilinear_form.apply2(intersection,
                                   *local_range_in_,
                                   *local_source_in_,
                                   *local_range_out_,
                                   *local_source_out_,
                                   bilinear_form_value_in_in_,
                                   bilinear_form_value_in_out_,
                                   bilinear_form_value_out_in_,
                                   bilinear_form_value_out_out_,
                                   param_);
        assert(bilinear_form_value_in_in_.rows() * bilinear_form_value_in_in_.cols() >= 1);
        assert(bilinear_form_value_in_out_.rows() * bilinear_form_value_in_out_.cols() >= 1);
        assert(bilinear_form_value_out_in_.rows() * bilinear_form_value_out_in_.cols() >= 1);
        assert(bilinear_form_value_out_out_.rows() * bilinear_form_value_out_out_.cols() >= 1);
        ret += bilinear_form_value_in_in_[0][0] + bilinear_form_value_in_out_[0][0] + bilinear_form_value_out_in_[0][0]
               + bilinear_form_value_out_out_[0][0];
      }
    }
    for (auto& data : bilinear_form_.intersection_data()) {
      const auto& local_bilinear_form = *std::get<0>(data);
      const auto& filter = *std::get<1>(data);
      auto apply_bilinear_form_to_one_side = [&](const auto& test_basis, const auto& ansatz_basis) {
        assert(test_basis.size(param_) == 1);
        assert(ansatz_basis.size(param_) == 1);
        // it does not matter which of the tmp matrices we pick here
        local_bilinear_form.apply2(intersection, test_basis, ansatz_basis, bilinear_form_value_in_in_, param_);
        assert(bilinear_form_value_in_in_.rows() * bilinear_form_value_in_in_.cols() >= 1);
        return bilinear_form_value_in_in_[0][0];
      };
      if (filter.contains(bilinear_form_.assembly_grid_view(), intersection)) {
        if (local_bilinear_form.inside())
          ret += apply_bilinear_form_to_one_side(*local_range_in_, *local_source_in_);
        else
          ret += apply_bilinear_form_to_one_side(*local_range_out_, *local_source_out_);
      }
    }
    return ret;
  } // ... compute_locally(...)

  const FieldType& result() const
  {
    return result_;
  }

  /// \}

  const BilinearFormType bilinear_form_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, r_r, r_rC, F>> range_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, s_r, s_rC, F>> source_;
  const XT::Common::Parameter param_;
  std::unique_ptr<typename RangeFunctionType::LocalFunctionType> local_range_in_;
  std::unique_ptr<typename RangeFunctionType::LocalFunctionType> local_range_out_;
  std::unique_ptr<typename SourceFunctionType::LocalFunctionType> local_source_in_;
  std::unique_ptr<typename SourceFunctionType::LocalFunctionType> local_source_out_;
  DynamicMatrix<FieldType> bilinear_form_value_in_in_;
  DynamicMatrix<FieldType> bilinear_form_value_in_out_;
  DynamicMatrix<FieldType> bilinear_form_value_out_in_;
  DynamicMatrix<FieldType> bilinear_form_value_out_out_;
  FieldType result_;
}; // class BilinearFormAssembler


template <class GridViewType,
          size_t s_r = 1,
          size_t r_r = s_r,
          size_t s_rC = 1,
          size_t r_rC = s_rC,
          class F = double,
          class SGV = GridViewType,
          class RGV = GridViewType>
auto make_bilinear_form(const GridViewType& grid_view,
                        const std::string& logging_prefix = "",
                        const std::array<bool, 3>& logging_state = {false, false, true})
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");
  return BilinearForm<GridViewType, s_r, s_rC, r_r, r_rC, F, SGV, RGV>(grid_view, logging_prefix, logging_state);
}


} // namespace GDT
} // namespace Dune

#include "matrix-based.hh"

#endif // DUNE_GDT_OPERATORS_BILINEAR_FORM_HH
