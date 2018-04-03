// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_BOUNDARY_HH
#define DUNE_GDT_OPERATORS_FV_BOUNDARY_HH

#include <dune/common/fvector.hh>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/boundaryinfo.hh>

namespace Dune {
namespace GDT {


template <class GridLayerImp, class RangeImp>
class LocalBoundaryValueInterface
{
public:
  using GridLayerType = GridLayerImp;
  using RangeType = RangeImp;
  using IntersectionType = typename GridLayerType::Intersection;
  using DomainFieldType = typename GridLayerType::ctype;
  static const size_t dimDomain = GridLayerType::dimension;
  using DomainType = FieldVector<DomainFieldType, dimDomain>;
  using BoundaryInfoType = XT::Grid::BoundaryInfo<IntersectionType>;
  using RangeFieldType = typename RangeType::field_type;

  LocalBoundaryValueInterface(const BoundaryInfoType& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual ~LocalBoundaryValueInterface() = default;

  virtual RangeType evaluate(const IntersectionType& intersection,
                             const DomainType& x,
                             const RangeType& u,
                             const XT::Common::Parameter& param = {}) const = 0;

protected:
  const BoundaryInfoType& boundary_info_;
};


template <class GridLayerImp, class RangeImp>
class LocalizableBoundaryValueInterface
{
public:
  using GridLayerType = GridLayerImp;
  using RangeType = RangeImp;
  using LocalBoundaryValueInterfaceType = LocalBoundaryValueInterface<GridLayerType, RangeType>;
  using BoundaryInfoType = typename LocalBoundaryValueInterfaceType::BoundaryInfoType;
  using DomainFieldType = typename LocalBoundaryValueInterfaceType::DomainFieldType;
  using DomainType = typename LocalBoundaryValueInterfaceType::DomainType;
  using EntityType = typename GridLayerType::template Codim<0>::Entity;
  using RangeFieldType = typename LocalBoundaryValueInterfaceType::RangeFieldType;

  LocalizableBoundaryValueInterface(const BoundaryInfoType& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual ~LocalizableBoundaryValueInterface() = default;

  virtual std::unique_ptr<LocalBoundaryValueInterfaceType> local_function(const EntityType& entity) const = 0;

protected:
  const BoundaryInfoType& boundary_info_;
};


template <class GridLayerType, class LocalizableFunctionType>
class LocalizableFunctionBasedLocalDirichletBoundaryValue
    : public LocalBoundaryValueInterface<GridLayerType, typename LocalizableFunctionType::RangeType>
{
  using BaseType = LocalBoundaryValueInterface<GridLayerType, typename LocalizableFunctionType::RangeType>;
  using LocalfunctionType = typename LocalizableFunctionType::LocalfunctionType;

public:
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::DomainType;
  using typename BaseType::IntersectionType;
  using typename BaseType::RangeType;

  LocalizableFunctionBasedLocalDirichletBoundaryValue(const BoundaryInfoType& boundary_info,
                                                      std::unique_ptr<LocalfunctionType>&& local_boundary_values)
    : BaseType(boundary_info)
    , local_boundary_values_(std::move(local_boundary_values))
  {
  }

  virtual RangeType
  evaluate(const IntersectionType& intersection, const DomainType& x, const RangeType& /*u*/) const override
  {
    if (boundary_info_.type(intersection) != XT::Grid::dirichlet_boundary)
      DUNE_THROW(Dune::NotImplemented, "This class can't handle boundary types other than Dirichlet!");
    return local_boundary_values_->evaluate(x);
  }

private:
  using BaseType::boundary_info_;
  const std::unique_ptr<LocalfunctionType> local_boundary_values_;
};


template <class GridLayerType, class LocalizableFunctionType>
class LocalizableFunctionBasedLocalizableDirichletBoundaryValue
    : public LocalizableBoundaryValueInterface<GridLayerType, typename LocalizableFunctionType::RangeType>
{
  using BaseType = LocalizableBoundaryValueInterface<GridLayerType, typename LocalizableFunctionType::RangeType>;
  using LocalBoundaryValueType =
      LocalizableFunctionBasedLocalDirichletBoundaryValue<GridLayerType, LocalizableFunctionType>;

public:
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalBoundaryValueInterfaceType;

  LocalizableFunctionBasedLocalizableDirichletBoundaryValue(const BoundaryInfoType& boundary_info,
                                                            const LocalizableFunctionType& boundary_values)
    : BaseType(boundary_info)
    , boundary_values_(boundary_values)
  {
  }

  virtual std::unique_ptr<LocalBoundaryValueInterfaceType> local_function(const EntityType& entity) const override
  {
    return LocalBoundaryValueType(boundary_info_, boundary_values_->local_function(entity));
  }

private:
  using BaseType::boundary_info_;
  const LocalizableFunctionType& boundary_values_;
};


template <class GridLayerType, class LocalizableFunctionType>
class LocalMomentModelBoundaryValue
    : public LocalBoundaryValueInterface<GridLayerType, typename LocalizableFunctionType::RangeType>
{
  using BaseType = LocalBoundaryValueInterface<GridLayerType, typename LocalizableFunctionType::RangeType>;
  using RangeFieldType = typename LocalizableFunctionType::RangeFieldType;

public:
  using MatrixType = typename Dune::XT::LA::CommonSparseOrDenseMatrixCsr<RangeFieldType>;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::DomainType;
  using typename BaseType::IntersectionType;
  using typename BaseType::RangeType;
  using BaseType::dimDomain;

  LocalMomentModelBoundaryValue(
      const BoundaryInfoType& boundary_info,
      const Dune::FieldVector<MatrixType, 2 * dimDomain>& reflection_matrices,
      std::unique_ptr<typename LocalizableFunctionType::LocalfunctionType>&& local_dirichlet_boundary_values)
    : BaseType(boundary_info)
    , reflection_matrices_(reflection_matrices)
    , local_dirichlet_boundary_values_(std::move(local_dirichlet_boundary_values))
  {
  }

  virtual RangeType evaluate(const IntersectionType& intersection,
                             const DomainType& x,
                             const RangeType& u,
                             const XT::Common::Parameter& /*param*/ = {}) const override
  {
    RangeType ret(0);
    const auto& boundary_type = boundary_info_.type(intersection);
    if (boundary_type == XT::Grid::absorbing_boundary)
      ret = u;
    else if (boundary_type == XT::Grid::dirichlet_boundary)
      ret = local_dirichlet_boundary_values_->evaluate(x);
    else if (boundary_type == XT::Grid::reflecting_boundary) {
      const auto direction = intersection.indexInInside();
      reflection_matrices_[direction].mv(u, ret);
    } else
      DUNE_THROW(Dune::NotImplemented,
                 "Boundary types other than Absorbing, Dirichlet or Reflecting are not implemented!");
    return ret;
  }

private:
  using BaseType::boundary_info_;
  const Dune::FieldVector<MatrixType, 2 * dimDomain>& reflection_matrices_;
  const std::unique_ptr<typename LocalizableFunctionType::LocalfunctionType> local_dirichlet_boundary_values_;
};


template <class GridLayerType, class BasisfunctionType, class LocalizableFunctionType>
class MomentModelBoundaryValue
    : public LocalizableBoundaryValueInterface<GridLayerType, typename BasisfunctionType::RangeType>
{
  static_assert(std::is_same<typename BasisfunctionType::RangeType, typename LocalizableFunctionType::RangeType>::value,
                "RangeTypes have to match!");
  using BaseType = LocalizableBoundaryValueInterface<GridLayerType, typename BasisfunctionType::RangeType>;

public:
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalBoundaryValueInterfaceType;
  using typename BaseType::RangeType;
  using LocalBoundaryValueType = LocalMomentModelBoundaryValue<GridLayerType, LocalizableFunctionType>;
  using MatrixType = typename LocalBoundaryValueType::MatrixType;
  static const size_t dimDomain = BasisfunctionType::dimDomain;
  static const size_t dimRange = BasisfunctionType::dimRange;


  MomentModelBoundaryValue(const BoundaryInfoType& boundary_info,
                           const BasisfunctionType& basis_functions,
                           const LocalizableFunctionType& dirichlet_boundary_values)
    : BaseType(boundary_info)
    , dirichlet_boundary_values_(dirichlet_boundary_values)
  {
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      DomainType n(0);
      n[ii] = -1.;
      reflection_matrices_[2 * ii] = MatrixType(basis_functions.reflection_matrix(n));
      n[ii] = 1.;
      reflection_matrices_[2 * ii + 1] = MatrixType(basis_functions.reflection_matrix(n));
    }
  }

  virtual std::unique_ptr<LocalBoundaryValueInterfaceType> local_function(const EntityType& entity) const override
  {
    return std::make_unique<LocalBoundaryValueType>(
        boundary_info_, reflection_matrices_, dirichlet_boundary_values_.local_function(entity));
  }

private:
  using BaseType::boundary_info_;
  Dune::FieldVector<MatrixType, 2 * dimDomain> reflection_matrices_;
  const LocalizableFunctionType& dirichlet_boundary_values_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_BOUNDARY_HH
