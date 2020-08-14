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
    , Logger(logging_prefix.empty() ? "gdt" : "gdt.bilinearform",
             logging_prefix.empty() ? "BilinearForm" : logging_prefix,
             /*logging_disabled=*/logging_prefix.empty())
    , grid_view_(grd_vw)
    , source_(src)
    , range_(rng)
    , local_source_(source_.local_function())
    , local_range_(range_.local_function())
    , bilinear_form_value_(1, 1, 0.)
    , result_(0)
  {
    LOG__(Logger, info) << Logger::logging_id << "(grid_view=" << &grd_vw << ", src=" << &src << ", rng=" << &rng << ")"
                        << std::endl;
  }

  BilinearForm(const ThisType& other)
    : BaseType(other)
    , Propagator(other)
    , Logger(other)
    , grid_view_(other.grid_view_)
    , source_(other.source_)
    , range_(other.range_)
    , local_source_(source_.local_function())
    , local_range_(range_.local_function())
    , bilinear_form_value_(other.bilinear_form_value_)
    , result_(other.result_)
  {
    for (auto& element_data : other.element_data_) {
      const auto& local_bilinear_form = *std::get<0>(element_data);
      auto& param = std::get<1>(element_data);
      const auto& filter = *std::get<2>(element_data);
      this->append(local_bilinear_form, param, filter);
    }
  } // BilinearForm(...)

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
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  ThisType&
  append(const LocalElementBilinearFormInterface<E, s_r, s_rC, SR, Result, r_r, r_rC, RR>& local_bilinear_form,
         const XT::Common::Parameter& param = {},
         const ElementFilterType& filter = XT::Grid::ApplyOn::AllElements<GV>())
  {
    LOG__(Logger, info) << Logger::logging_id << ".append(local_bilinear_form=" << &local_bilinear_form
                        << ", param=" << param << ", filter=" << &filter << ")" << std::endl;
    element_data_.emplace_back(local_bilinear_form.copy(), param, std::unique_ptr<ElementFilterType>(filter.copy()));
    return *this;
  }

  ThisType&
  operator+=(const LocalElementBilinearFormInterface<E, s_r, s_rC, SR, Result, r_r, r_rC, RR>& local_bilinear_form)
  {
    this->append(local_bilinear_form);
    return *this;
  }

  ThisType& operator+=(std::tuple<const LocalElementBilinearFormInterface<E, s_r, s_rC, SR, Result, r_r, r_rC, RR>&,
                                  const XT::Common::Parameter&> local_bilinear_form__param)
  {
    this->append(std::get<0>(local_bilinear_form__param), std::get<1>(local_bilinear_form__param));
    return *this;
  }

  ThisType& operator+=(std::tuple<const LocalElementBilinearFormInterface<E, s_r, s_rC, SR, Result, r_r, r_rC, RR>&,
                                  const ElementFilterType&> local_bilinear_form__filter)
  {
    this->append(std::get<0>(local_bilinear_form__filter), {}, std::get<1>(local_bilinear_form__filter));
    return *this;
  }

  ThisType& operator+=(std::tuple<const LocalElementBilinearFormInterface<E, s_r, s_rC, SR, Result, r_r, r_rC, RR>&,
                                  const XT::Common::Parameter&,
                                  const ElementFilterType&> local_bilinear_form__param__filter)
  {
    this->append(std::get<0>(local_bilinear_form__param__filter),
                 std::get<1>(local_bilinear_form__param__filter),
                 std::get<2>(local_bilinear_form__param__filter));
    return *this;
  }

  void prepare() override final
  {
    result_ = 0;
  }

  Result compute_locally(const E& element)
  {
    Result ret = 0;
    local_source_->bind(element);
    local_range_->bind(element);
    DUNE_THROW_IF(local_source_->size() != 1,
                  Exceptions::bilinear_form_error,
                  "local_source_->size() = " << local_source_->size());
    DUNE_THROW_IF(
        local_range_->size() != 1, Exceptions::bilinear_form_error, "local_range_->size() = " << local_range_->size());
    for (auto& element_data : element_data_) {
      const auto& local_bilinear_form = *std::get<0>(element_data);
      auto& param = std::get<1>(element_data);
      const auto& filter = *std::get<2>(element_data);
      if (filter.contains(this->grid_view_, element)) {
        local_bilinear_form.apply2(*local_source_, *local_range_, bilinear_form_value_, param);
        assert(bilinear_form_value_.rows() >= 1);
        assert(bilinear_form_value_.cols() >= 1);
        ret += bilinear_form_value_[0][0];
      }
    }
    return ret;
  }

  void apply_local(const E& element) override final
  {
    result_ += compute_locally(element);
  }

  void apply_local(const I& /*intersection*/, const E& /*inside_element*/, const E& /*outside_element*/) override final
  {}

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
  XT::Functions::GridFunction<E, s_r, s_rC, SR> source_;
  XT::Functions::GridFunction<E, r_r, r_rC, RR> range_;
  std::unique_ptr<typename SourceType::LocalFunctionType> local_source_;
  std::unique_ptr<typename RangeType::LocalFunctionType> local_range_;
  DynamicMatrix<Result> bilinear_form_value_;
  ResultType result_;
  std::list<std::tuple<std::unique_ptr<LocalElementBilinearFormInterface<E, s_r, s_rC, SR, ResultType, r_r, r_rC, RR>>,
                       XT::Common::Parameter,
                       std::unique_ptr<ElementFilterType>>>
      element_data_;
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
