// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_TOOLS_DISCRETEVALUED_GRID_FUNCTION_HH
#define DUNE_GDT_TOOLS_DISCRETEVALUED_GRID_FUNCTION_HH

#include <memory>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/vector_less.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Wrapper for the map of reconstructed values that fulfills the XT::Functions::LocalizableFunctionInterface
 */
template <class GV, size_t rangeDim, size_t rangeDimCols, class RangeField>
class DiscreteValuedGridFunction
  : public XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, rangeDim, rangeDimCols, RangeField>
{
  using BaseType =
      XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, rangeDim, rangeDimCols, RangeField>;
  using ThisType = DiscreteValuedGridFunction;

public:
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;

  using IndexSetType = typename GV::IndexSet;
  using E = XT::Grid::extract_entity_t<GV>;
  using typename BaseType::LocalFunctionType;
  using typename BaseType::R;

private:
  class DiscreteValuedLocalFunction : public LocalFunctionType
  {
    using BaseType = LocalFunctionType;

  public:
    using typename BaseType::DomainType;
    using typename BaseType::DynamicRangeType;
    using typename BaseType::E;
    using typename BaseType::RangeReturnType;
    using LocalFunctionValuesType = std::map<DomainType, RangeReturnType, XT::Common::FieldVectorLess>;

    DiscreteValuedLocalFunction(std::vector<LocalFunctionValuesType>& values, const IndexSetType& index_set)
      : values_(values)
      , index_set_(index_set)
    {}

    int order(const XT::Common::Parameter& /*mu*/ = {}) const override
    {
      DUNE_THROW(Dune::InvalidStateException, "This function can't be integrated!");
      return 2;
    }

    using BaseType::element;

    void post_bind(const E& elem) override final
    {
      local_values_ = &(values_[index_set_.index(elem)]);
    }

    RangeReturnType evaluate(const DomainType& xx, const XT::Common::Parameter& /*param*/) const override final
    {
      try {
        return local_values_->at(xx);
      } catch (const std::out_of_range& /*e*/) {
        DUNE_THROW(Dune::RangeError,
                   "There are no values for local coord "
                       << XT::Common::to_string(xx) << " (global coord "
                       << XT::Common::to_string(element().geometry().global(xx)) << ") on entity "
                       << XT::Common::to_string(element().geometry().center()) << " in this function!");
      }
      return RangeReturnType{};
    }

    void
    evaluate(const DomainType& xx, DynamicRangeType& ret, const XT::Common::Parameter& /*param*/) const override final
    {
      try {
        const auto& vals = local_values_->at(xx);
        for (size_t ii = 0; ii < r; ++ii)
          ret[ii] = vals[ii];
      } catch (const std::out_of_range& /*e*/) {
        DUNE_THROW(Dune::RangeError,
                   "There are no values for local coord "
                       << XT::Common::to_string(xx) << " (global coord "
                       << XT::Common::to_string(element().geometry().global(xx)) << ") on entity "
                       << XT::Common::to_string(element().geometry().center()) << " in this function!");
      }
    }

  private:
    std::vector<LocalFunctionValuesType>& values_;
    const IndexSetType& index_set_;
    LocalFunctionValuesType* local_values_;
  };

public:
  using LocalFunctionValuesType = typename DiscreteValuedLocalFunction::LocalFunctionValuesType;

  static const bool available = true;

  DiscreteValuedGridFunction(const GV& grid_view, std::vector<LocalFunctionValuesType>& values)
    : index_set_(grid_view.indexSet())
    , values_(values)
  {
    assert(grid_view.size(0) >= 0);
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<DiscreteValuedLocalFunction>(values_, index_set_);
  }

  LocalFunctionValuesType& local_values(const E& elem)
  {
    return values_[index_set_.index(elem)];
  }

private:
  const IndexSetType& index_set_;
  std::vector<LocalFunctionValuesType>& values_;
}; // class DiscreteValuedGridFunction


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_DISCRETEVALUED_GRID_FUNCTION_HH
