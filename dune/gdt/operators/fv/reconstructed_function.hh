// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTED_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTED_HH

#include <memory>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/functions/interfaces/localizable-function.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Wrapper for the map of reconstructed values that fulfills the XT::Functions::LocalizableFunctionInterface
 */
template <class GridViewImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class ReconstructedLocalizableFunction
    : public XT::Functions::LocalizableFunctionInterface<typename GridViewImp::template Codim<0>::Entity,
                                                         DomainFieldImp,
                                                         domainDim,
                                                         RangeFieldImp,
                                                         rangeDim,
                                                         rangeDimCols>
{
  typedef XT::Functions::LocalizableFunctionInterface<typename GridViewImp::template Codim<0>::Entity,
                                                      DomainFieldImp,
                                                      domainDim,
                                                      RangeFieldImp,
                                                      rangeDim,
                                                      rangeDimCols>
      BaseType;

public:
  static const constexpr size_t dimDomain = BaseType::dimDomain;
  static const constexpr size_t dimRange = BaseType::dimRange;
  static const constexpr size_t dimRangeCols = BaseType::dimRangeCols;

  typedef GridViewImp GridLayerType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
  typedef typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate LocalCoordinateType;

private:
  class ReconstructedLocalfunction
      : public XT::Functions::
            LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
  {
    typedef typename XT::Functions::
        LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
            BaseType;

  public:
    ReconstructedLocalfunction(const EntityType& ent,
                               const std::map<LocalCoordinateType, RangeType, XT::Common::FieldVectorLess>& values)
      : BaseType(ent)
      , values_(values)
    {
    }

    virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const
    {
      DUNE_THROW(Dune::InvalidStateException, "This function can't be integrated!");
      return 2;
    }

    using BaseType::entity;

    virtual void evaluate(const DomainType& xx, RangeType& ret, const XT::Common::Parameter& /*param*/) const
    {
      try {
        ret = values_.at(xx);
      } catch (const std::out_of_range& /*e*/) {
        DUNE_THROW(Dune::RangeError,
                   "There are no values for local coord " << XT::Common::to_string(xx) << " (global coord "
                                                          << XT::Common::to_string(entity().geometry().global(xx))
                                                          << ") on entity "
                                                          << XT::Common::to_string(entity().geometry().center())
                                                          << " in this function!");
      }
    }

    virtual void
    jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/, const XT::Common::Parameter& /*param*/) const
    {
      DUNE_THROW(Dune::NotImplemented, "");
    }

  private:
    const std::map<LocalCoordinateType, RangeType, XT::Common::FieldVectorLess>& values_;
  };

public:
  typedef ReconstructedLocalfunction LocalfunctionType;
  typedef XT::Functions::
      LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
          LocalfunctionInterfaceType;

  static const bool available = true;

  ReconstructedLocalizableFunction(
      const GridLayerType& grid_layer,
      const std::vector<std::map<LocalCoordinateType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values =
          {})
    : index_set_(grid_layer.indexSet())
    , reconstructed_values_(reconstructed_values)
  {
  }

  virtual std::unique_ptr<LocalfunctionInterfaceType> local_function(const EntityType& entity) const
  {
    return std::make_unique<LocalfunctionType>(entity, reconstructed_values_[index_set_.index(entity)]);
  }

  static std::string static_id()
  {
    return "reconstructed localizable function";
  }

  virtual std::string type() const
  {
    return "reconstructed localizable function";
  }

  virtual std::string name() const
  {
    return "reconstructed localizable function";
  }

  void resize(const size_t size)
  {
    reconstructed_values_.resize(size);
  }

private:
  const IndexSetType& index_set_;
  const std::vector<std::map<LocalCoordinateType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values_;
}; // class ReconstructedLocalizableFunction


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTED_HH
