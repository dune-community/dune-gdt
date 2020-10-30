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

#include <dune/gdt/local/assembler/bilinear-form-accumulators.hh>
#include <dune/gdt/local/bilinear-forms/interfaces.hh>

namespace Dune {
namespace GDT {


template <class GridView,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceRangeField = double,
          class Result = double,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeRangeField = SourceRangeField>
class BilinearForm
  : public XT::Grid::ElementAndIntersectionFunctor<GridView>
  , XT::Common::ThreadResultPropagator<BilinearForm<GridView,
                                                    source_range_dim,
                                                    source_range_dim_cols,
                                                    SourceRangeField,
                                                    Result,
                                                    range_range_dim,
                                                    range_range_dim_cols,
                                                    RangeRangeField>,
                                       Result>
  , public XT::Common::WithLogger<BilinearForm<GridView,
                                               source_range_dim,
                                               source_range_dim_cols,
                                               SourceRangeField,
                                               Result,
                                               range_range_dim,
                                               range_range_dim_cols,
                                               RangeRangeField>>
{
  static_assert(XT::Grid::is_view<GridView>::value, "");

  using ThisType = BilinearForm;
  using BaseType = XT::Grid::ElementAndIntersectionFunctor<GridView>;
  using Propagator = XT::Common::ThreadResultPropagator<BilinearForm<GridView,
                                                                     source_range_dim,
                                                                     source_range_dim_cols,
                                                                     SourceRangeField,
                                                                     Result,
                                                                     range_range_dim,
                                                                     range_range_dim_cols,
                                                                     RangeRangeField>,
                                                        Result>;
  friend Propagator;

  using Logger = XT::Common::WithLogger<BilinearForm<GridView,
                                                     source_range_dim,
                                                     source_range_dim_cols,
                                                     SourceRangeField,
                                                     Result,
                                                     range_range_dim,
                                                     range_range_dim_cols,
                                                     RangeRangeField>>;

public:
  using GV = GridView;
  using typename BaseType::E;
  using typename BaseType::GridViewType;
  using typename BaseType::I;
  using ResultType = Result;

  static const constexpr size_t s_r = source_range_dim;
  static const constexpr size_t s_rC = source_range_dim_cols;
  using SR = SourceRangeField;
  using SourceType = XT::Functions::GridFunctionInterface<E, s_r, s_rC, SR>;

  static const constexpr size_t r_r = range_range_dim;
  static const constexpr size_t r_rC = range_range_dim_cols;
  using RR = RangeRangeField;
  using RangeType = XT::Functions::GridFunctionInterface<E, r_r, r_rC, RR>;

  using ElementFilterType = XT::Grid::ElementFilter<GV>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<GV>;

  BilinearForm(GridViewType grd_vw,
               XT::Functions::GridFunction<E, s_r, s_rC, SR> src,
               XT::Functions::GridFunction<E, r_r, r_rC, RR> rng,
               const std::string& logging_prefix = "")
    : BaseType()
    , Propagator(this)
    , Logger(logging_prefix.empty() ? "BilinearForm" : logging_prefix,
             /*logging_disabled=*/logging_prefix.empty())
    , grid_view_(grd_vw)
    , source_(src.copy_as_grid_function())
    , range_(rng.copy_as_grid_function())
    , local_source_in_(source_->local_function())
    , local_source_out_(source_->local_function())
    , local_range_in_(range_->local_function())
    , local_range_out_(range_->local_function())
    , bilinear_form_value_in_in_(1, 1, 0.)
    , bilinear_form_value_in_out_(1, 1, 0.)
    , bilinear_form_value_out_in_(1, 1, 0.)
    , bilinear_form_value_out_out_(1, 1, 0.)
    , result_(0)
  {
    LOG__(Logger, info) << Logger::logger.prefix << "BilinearForm(grid_view=" << &grd_vw << ", src=" << &src
                        << ", rng=" << &rng << ")" << std::endl;
  }

  BilinearForm(const ThisType& other)
    : BaseType(other)
    , Propagator(other)
    , Logger(other)
    , grid_view_(other.grid_view_)
    , source_(other.source_->copy_as_grid_function())
    , range_(other.range_->copy_as_grid_function())
    , local_source_in_(source_->local_function())
    , local_source_out_(source_->local_function())
    , local_range_in_(range_->local_function())
    , local_range_out_(range_->local_function())
    , bilinear_form_value_in_in_(other.bilinear_form_value_in_in_)
    , bilinear_form_value_in_out_(other.bilinear_form_value_in_out_)
    , bilinear_form_value_out_in_(other.bilinear_form_value_out_in_)
    , bilinear_form_value_out_out_(other.bilinear_form_value_out_out_)
    , result_(other.result_)
  {
    for (auto& data : other.element_data_) {
      const auto& local_bilinear_form = *std::get<0>(data);
      auto& param = std::get<1>(data);
      const auto& filter = *std::get<2>(data);
      this->append(local_bilinear_form, param, filter);
    }
    for (auto& data : other.coupling_intersection_data_) {
      const auto& local_bilinear_form = *std::get<0>(data);
      auto& param = std::get<1>(data);
      const auto& filter = *std::get<2>(data);
      this->append(local_bilinear_form, param, filter);
    }
    for (auto& data : other.intersection_data_) {
      const auto& local_bilinear_form = *std::get<0>(data);
      auto& param = std::get<1>(data);
      const auto& filter = *std::get<2>(data);
      this->append(local_bilinear_form, param, filter);
    }
  } // BilinearForm(...)

  BilinearForm(ThisType&&) = default;

  BaseType* copy() override final
  {
    return Propagator::copy_imp();
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const SourceType& source() const
  {
    return *source_;
  }

  const RangeType& range() const
  {
    return *range_;
  }

  /// \name These methods allow to append local element bilinear forms
  ///\{

  ThisType&
  append(const LocalElementBilinearFormInterface<E, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>& local_bilinear_form,
         const XT::Common::Parameter& param = {},
         const ElementFilterType& filter = XT::Grid::ApplyOn::AllElements<GV>())
  {
    LOG__(Logger, info) << Logger::logger.prefix << ".append(local_element_bilinear_form=" << &local_bilinear_form
                        << ", param=" << param << ", filter=" << &filter << ")" << std::endl;
    element_data_.emplace_back(local_bilinear_form.copy(), param, std::unique_ptr<ElementFilterType>(filter.copy()));
    return *this;
  }

  ThisType&
  operator+=(const LocalElementBilinearFormInterface<E, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>& local_bilinear_form)
  {
    this->append(local_bilinear_form);
    return *this;
  }

  ThisType&
  operator+=(const std::tuple<const LocalElementBilinearFormInterface<E, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
                              const XT::Common::Parameter&>& local_bilinear_form__param)
  {
    this->append(std::get<0>(local_bilinear_form__param), std::get<1>(local_bilinear_form__param));
    return *this;
  }

  ThisType&
  operator+=(const std::tuple<const LocalElementBilinearFormInterface<E, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
                              const ElementFilterType&>& local_bilinear_form__filter)
  {
    this->append(std::get<0>(local_bilinear_form__filter), {}, std::get<1>(local_bilinear_form__filter));
    return *this;
  }

  ThisType&
  operator+=(const std::tuple<const LocalElementBilinearFormInterface<E, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
                              const XT::Common::Parameter&,
                              const ElementFilterType&>& local_bilinear_form__param__filter)
  {
    this->append(std::get<0>(local_bilinear_form__param__filter),
                 std::get<1>(local_bilinear_form__param__filter),
                 std::get<2>(local_bilinear_form__param__filter));
    return *this;
  }

  ///\}
  /// \name These methods allow to append local coupling intersection bilinear forms
  ///\{

  ThisType& append(const LocalCouplingIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&
                       local_bilinear_form,
                   const XT::Common::Parameter& param = {},
                   const IntersectionFilterType& filter = XT::Grid::ApplyOn::AllIntersections<GV>())
  {
    LOG__(Logger, info) << Logger::logger.prefix
                        << ".append(local_coupling_intersection_bilinear_form=" << &local_bilinear_form
                        << ", param=" << param << ", filter=" << &filter << ")" << std::endl;
    coupling_intersection_data_.emplace_back(
        local_bilinear_form.copy(), param, std::unique_ptr<IntersectionFilterType>(filter.copy()));
    return *this;
  }

  ThisType&
  operator+=(const LocalCouplingIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&
                 local_bilinear_form)
  {
    this->append(local_bilinear_form);
    return *this;
  }

  ThisType&
  operator+=(const std::tuple<
             const LocalCouplingIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
             const XT::Common::Parameter&>& local_bilinear_form__param)
  {
    this->append(std::get<0>(local_bilinear_form__param), std::get<1>(local_bilinear_form__param));
    return *this;
  }

  ThisType&
  operator+=(const std::tuple<
             const LocalCouplingIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
             const IntersectionFilterType&>& local_bilinear_form__filter)
  {
    this->append(std::get<0>(local_bilinear_form__filter), {}, std::get<1>(local_bilinear_form__filter));
    return *this;
  }

  ThisType&
  operator+=(const std::tuple<
             const LocalCouplingIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
             const XT::Common::Parameter&,
             const IntersectionFilterType&>& local_bilinear_form__param__filter)
  {
    this->append(std::get<0>(local_bilinear_form__param__filter),
                 std::get<1>(local_bilinear_form__param__filter),
                 std::get<2>(local_bilinear_form__param__filter));
    return *this;
  }

  ///\}
  /// \name These methods allow to append local intersection bilinear forms
  ///\{

  ThisType&
  append(const LocalIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>& local_bilinear_form,
         const XT::Common::Parameter& param = {},
         const IntersectionFilterType& filter = XT::Grid::ApplyOn::AllIntersections<GV>())
  {
    LOG__(Logger, info) << Logger::logger.prefix << ".append(local_intersection_bilinear_form=" << &local_bilinear_form
                        << ", param=" << param << ", filter=" << &filter << ")" << std::endl;
    intersection_data_.emplace_back(
        local_bilinear_form.copy(), param, std::unique_ptr<IntersectionFilterType>(filter.copy()));
    return *this;
  }

  ThisType& operator+=(
      const LocalIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>& local_bilinear_form)
  {
    this->append(local_bilinear_form);
    return *this;
  }

  ThisType& operator+=(
      const std::tuple<const LocalIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
                       const XT::Common::Parameter&>& local_bilinear_form__param)
  {
    this->append(std::get<0>(local_bilinear_form__param), std::get<1>(local_bilinear_form__param));
    return *this;
  }

  ThisType& operator+=(
      const std::tuple<const LocalIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
                       const IntersectionFilterType&>& local_bilinear_form__filter)
  {
    this->append(std::get<0>(local_bilinear_form__filter), {}, std::get<1>(local_bilinear_form__filter));
    return *this;
  }

  ThisType& operator+=(
      const std::tuple<const LocalIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>&,
                       const XT::Common::Parameter&,
                       const IntersectionFilterType&>& local_bilinear_form__param__filter)
  {
    this->append(std::get<0>(local_bilinear_form__param__filter),
                 std::get<1>(local_bilinear_form__param__filter),
                 std::get<2>(local_bilinear_form__param__filter));
    return *this;
  }

  ///\}

  void prepare() override final
  {
    result_ = 0;
  }

  /// \brief Variant of compute_locally to apply all local element bilinear forms
  ResultType compute_locally(const E& element)
  {
    Result ret = 0;
    local_source_in_->bind(element);
    local_range_in_->bind(element);
    for (auto& data : element_data_) {
      const auto& local_bilinear_form = *std::get<0>(data);
      auto& param = std::get<1>(data);
      const auto& filter = *std::get<2>(data);
      if (filter.contains(this->grid_view_, element)) {
        assert(local_source_in_->size(param) == 1);
        assert(local_range_in_->size(param) == 1);
        local_bilinear_form.apply2(*local_range_in_, *local_source_in_, bilinear_form_value_in_in_, param);
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
  ResultType compute_locally(const I& intersection, const E& inside_element, const E& outside_element)
  {
    Result ret = 0;
    local_source_in_->bind(inside_element);
    local_range_in_->bind(inside_element);
    local_source_out_->bind(outside_element);
    local_range_out_->bind(outside_element);
    for (auto& data : coupling_intersection_data_) {
      const auto& local_bilinear_form = *std::get<0>(data);
      auto& param = std::get<1>(data);
      const auto& filter = *std::get<2>(data);
      if (filter.contains(this->grid_view_, intersection)) {
        assert(local_source_in_->size(param) == 1);
        assert(local_source_out_->size(param) == 1);
        assert(local_range_in_->size(param) == 1);
        assert(local_range_out_->size(param) == 1);
        local_bilinear_form.apply2(intersection,
                                   *local_range_in_,
                                   *local_source_in_,
                                   *local_range_out_,
                                   *local_source_out_,
                                   bilinear_form_value_in_in_,
                                   bilinear_form_value_in_out_,
                                   bilinear_form_value_out_in_,
                                   bilinear_form_value_out_out_,
                                   param);
        assert(bilinear_form_value_in_in_.rows() * bilinear_form_value_in_in_.cols() >= 1);
        assert(bilinear_form_value_in_out_.rows() * bilinear_form_value_in_out_.cols() >= 1);
        assert(bilinear_form_value_out_in_.rows() * bilinear_form_value_out_in_.cols() >= 1);
        assert(bilinear_form_value_out_out_.rows() * bilinear_form_value_out_out_.cols() >= 1);
        ret += bilinear_form_value_in_in_[0][0] + bilinear_form_value_in_out_[0][0] + bilinear_form_value_out_in_[0][0]
               + bilinear_form_value_out_out_[0][0];
      }
    }
    for (auto& data : intersection_data_) {
      const auto& local_bilinear_form = *std::get<0>(data);
      auto& param = std::get<1>(data);
      const auto& filter = *std::get<2>(data);
      auto apply_bilinear_form_to_one_side = [&](const auto& test_basis, const auto& ansatz_basis) {
        assert(test_basis.size(param) == 1);
        assert(ansatz_basis.size(param) == 1);
        // it does not matter which of the tmp matrices we pick here
        local_bilinear_form.apply2(intersection, test_basis, ansatz_basis, bilinear_form_value_in_in_, param);
        assert(bilinear_form_value_in_in_.rows() * bilinear_form_value_in_in_.cols() >= 1);
        return bilinear_form_value_in_in_[0][0];
      };
      if (filter.contains(this->grid_view_, intersection)) {
        if (local_bilinear_form.inside())
          ret += apply_bilinear_form_to_one_side(*local_range_in_, *local_source_in_);
        else
          ret += apply_bilinear_form_to_one_side(*local_range_out_, *local_source_out_);
      }
    }
    return ret;
  } // ... compute_locally(...)

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

  const ResultType& result() const
  {
    return result_;
  }

  Result apply2(const bool use_tbb = false)
  {
    this->prepare();
    XT::Grid::Walker<GV> walker(this->grid_view_);
    walker.append(*this);
    walker.walk(use_tbb);
    return this->result_;
  }

protected:
  void set_result(Result res)
  {
    result_ = res;
  }

  const GridViewType grid_view_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, s_r, s_rC, SR>> source_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, r_r, r_rC, RR>> range_;
  std::unique_ptr<typename SourceType::LocalFunctionType> local_source_in_;
  std::unique_ptr<typename SourceType::LocalFunctionType> local_source_out_;
  std::unique_ptr<typename RangeType::LocalFunctionType> local_range_in_;
  std::unique_ptr<typename RangeType::LocalFunctionType> local_range_out_;
  DynamicMatrix<Result> bilinear_form_value_in_in_;
  DynamicMatrix<Result> bilinear_form_value_in_out_;
  DynamicMatrix<Result> bilinear_form_value_out_in_;
  DynamicMatrix<Result> bilinear_form_value_out_out_;
  ResultType result_;
  std::list<std::tuple<std::unique_ptr<LocalElementBilinearFormInterface<E, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>>,
                       XT::Common::Parameter,
                       std::unique_ptr<ElementFilterType>>>
      element_data_;
  std::list<std::tuple<
      std::unique_ptr<LocalCouplingIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>>,
      XT::Common::Parameter,
      std::unique_ptr<IntersectionFilterType>>>
      coupling_intersection_data_;
  std::list<
      std::tuple<std::unique_ptr<LocalIntersectionBilinearFormInterface<I, r_r, r_rC, RR, ResultType, s_r, s_rC, SR>>,
                 XT::Common::Parameter,
                 std::unique_ptr<IntersectionFilterType>>>
      intersection_data_;
}; // class BilinearForm


template <class GV, size_t s_r, size_t s_rC, class SR, size_t r_r, size_t r_rC, class RR>
std::enable_if_t<XT::Grid::is_view<GV>::value, BilinearForm<GV, s_r, s_rC, SR, double, r_r, r_rC, RR>>
make_bilinear_form(GV grid_view,
                   const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, s_r, s_rC, SR>& source,
                   const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r_r, r_rC, RR>& range,
                   const std::string& logging_prefix = "")
{
  return BilinearForm<GV, s_r, s_rC, SR, double, r_r, r_rC, RR>(grid_view, source, range, logging_prefix);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_BILINEAR_FORM_HH
