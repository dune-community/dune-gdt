// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2011 - 2018)
//   Kirsten Weber   (2013)
//   René Fritze     (2014 - 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/la/container.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/local/discretefunction.hh>
#include <dune/gdt/discretefunction/dof-vector.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace internal {


template <class Vector, class GridView>
struct AssertArgumentsOfConstDiscreteFunction
{
  static_assert(XT::LA::is_vector<Vector>::value, "");
  static_assert(XT::Grid::is_view<GridView>::value, "");
  using E = XT::Grid::extract_entity_t<GridView>;
};


} // namespace internal


template <class Vector, class GridView, size_t range_dim = 1, size_t range_dim_cols = 1, class RangeField = double>
class ConstDiscreteFunction
  : public XT::Functions::GridFunctionInterface<
        typename internal::AssertArgumentsOfConstDiscreteFunction<Vector, GridView>::E,
        range_dim,
        range_dim_cols,
        RangeField>
{
  using BaseType = XT::Functions::GridFunctionInterface<
      typename internal::AssertArgumentsOfConstDiscreteFunction<Vector, GridView>::E,
      range_dim,
      range_dim_cols,
      RangeField>;
  using ThisType = ConstDiscreteFunction;

public:
  using ConstDofVectorType = ConstDofVector<Vector, GridView>;
  using SpaceType = SpaceInterface<GridView, range_dim, range_dim_cols, RangeField>;
  using VectorType = Vector;
  using ConstLocalDiscreteFunctionType =
      ConstLocalDiscreteFunction<Vector, GridView, range_dim, range_dim_cols, RangeField>;

  using typename BaseType::ElementType;
  using typename BaseType::LocalFunctionType;

  ConstDiscreteFunction(const SpaceType& spc,
                        const VectorType& vector,
                        const std::string nm = "dune.gdt.constdiscretefunction")
    : space_(spc)
    , dofs_(space_.mapper(), vector)
    , name_(nm)
  {}

  ConstDiscreteFunction(const ThisType&) = default;
  ConstDiscreteFunction(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const SpaceType& space() const
  {
    return space_;
  }

  const ConstDofVectorType& dofs() const
  {
    return dofs_;
  }

  std::unique_ptr<ConstLocalDiscreteFunctionType> local_discrete_function() const
  {
    return std::make_unique<ConstLocalDiscreteFunctionType>(space_, dofs_);
  }

  std::unique_ptr<ConstLocalDiscreteFunctionType> local_discrete_function(const ElementType& grid_element) const
  {
    auto ldf = local_discrete_function();
    ldf->bind(grid_element);
    return ldf;
  }

  /**
   * \name ``These methods are required by XT::Functions::GridFunctionInterface.''
   * \{
   */

  std::string name() const override final
  {
    return name_;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return local_discrete_function();
  }

  /**
   * \}
   */

  using BaseType::visualize;
  using BaseType::visualize_gradient;

  /**
   * \brief Visualizes the function using Dune::XT::Functions::GridFunctionInterface::visualize on the grid view
   *        associated with the space.
   * \sa    Dune::XT::Functions::GridFunctionInterface::visualize
   * \note  Subsampling is enabled by default for functions of order greater than one.
   */
  void visualize(const std::string filename,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw,
                 const XT::Common::Parameter& param = {}) const
  {
    const bool subsampling =
        param.has_key("subsampling") ? static_cast<bool>(param.get("subsampling")[0]) : (space_.max_polorder() > 1);
    this->visualize(space_.grid_view(), filename, subsampling, vtk_output_type, param);
  }

  /**
   * \brief Visualizes the function using Dune::XT::Functions::GridFunctionInterface::visualize on the grid view
   *        associated with the space.
   * \sa    Dune::XT::Functions::GridFunctionInterface::visualize
   * \note  Subsampling is enabled by default for functions of order greater than one.
   */
  void visualize_gradient(const std::string filename,
                          const VTK::OutputType vtk_output_type = VTK::appendedraw,
                          const XT::Common::Parameter& param = {}) const
  {
    const bool subsampling =
        param.has_key("subsampling") ? static_cast<bool>(param.get("subsampling")[0]) : (space_.max_polorder() > 1);
    this->visualize_gradient(space_.grid_view(), filename, subsampling, vtk_output_type, param);
  }


protected:
  const SpaceType& space_;

private:
  const ConstDofVectorType dofs_;
  const std::string name_;
}; // class ConstDiscreteFunction


template <class V, class GV, size_t r, size_t rC, class R>
ConstDiscreteFunction<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, R>
make_discrete_function(const SpaceInterface<GV, r, rC, R>& space,
                       const XT::LA::VectorInterface<V>& vector,
                       const std::string nm = "dune.gdt.constdiscretefunction")
{
  return ConstDiscreteFunction<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, R>(
      space, vector.as_imp(), nm);
}


template <class Vector, class GridView, size_t range_dim = 1, size_t range_dim_cols = 1, class RangeField = double>
class DiscreteFunction
  : XT::Common::StorageProvider<Vector>
  , public ConstDiscreteFunction<Vector, GridView, range_dim, range_dim_cols, RangeField>
{
  using ThisType = DiscreteFunction;
  using VectorStorage = XT::Common::StorageProvider<Vector>;
  using BaseType = ConstDiscreteFunction<Vector, GridView, range_dim, range_dim_cols, RangeField>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::SpaceType;
  using typename BaseType::VectorType;

  using DofVectorType = DofVector<Vector, GridView>;
  using LocalDiscreteFunctionType = LocalDiscreteFunction<Vector, GridView, range_dim, range_dim_cols, RangeField>;

  DiscreteFunction(const SpaceType& spc, VectorType& vector, const std::string nm = "dune.gdt.discretefunction")
    : VectorStorage(vector)
    , BaseType(spc, VectorStorage::access(), nm)
    , dofs_(space_.mapper(), VectorStorage::access())
  {}

  DiscreteFunction(const SpaceType& spc, VectorType&& vector, const std::string nm = "dune.gdt.discretefunction")
    : VectorStorage(new VectorType(std::move(vector)))
    , BaseType(spc, VectorStorage::access(), nm)
    , dofs_(space_.mapper(), VectorStorage::access())
  {}

  DiscreteFunction(const SpaceType& spc, const std::string nm = "dune.gdt.discretefunction")
    : VectorStorage(new VectorType(spc.mapper().size(), 0.))
    , BaseType(spc, VectorStorage::access(), nm)
    , dofs_(space_.mapper(), VectorStorage::access())
  {}

  DiscreteFunction(const ThisType& other)
    : VectorStorage(new VectorType(other.access()))
    , BaseType(other.space(), VectorStorage::access(), other.name())
    , dofs_(other.space().mapper(), VectorStorage::access())
  {}

  DiscreteFunction(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  using BaseType::dofs;

  DofVectorType& dofs()
  {
    return dofs_;
  }

  using BaseType::local_discrete_function;

  std::unique_ptr<LocalDiscreteFunctionType> local_discrete_function()
  {
    return std::make_unique<LocalDiscreteFunctionType>(space_, dofs_);
  }

  std::unique_ptr<LocalDiscreteFunctionType> local_discrete_function(const ElementType& grid_element)
  {
    auto ldf = local_discrete_function();
    ldf->bind(grid_element);
    return ldf;
  }
  /**
   * \}
   */

  ThisType& operator+=(const BaseType& other)
  {
    dofs().vector() += other.dofs().vector();
    return *this;
  }

  ThisType& operator-=(const BaseType& other)
  {
    dofs().vector() -= other.dofs().vector();
    return *this;
  }

private:
  using BaseType::space_;
  DofVectorType dofs_;
}; // class DiscreteFunction


template <class V, class GV, size_t r, size_t rC, class R>
DiscreteFunction<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, R>
make_discrete_function(const SpaceInterface<GV, r, rC, R>& space,
                       XT::LA::VectorInterface<V>& vector,
                       const std::string nm = "dune.gdt.discretefunction")
{
  return DiscreteFunction<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, R>(space, vector.as_imp(), nm);
}


template <class V, class GV, size_t r, size_t rC, class R>
DiscreteFunction<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, R>
make_discrete_function(const SpaceInterface<GV, r, rC, R>& space,
                       XT::LA::VectorInterface<V>&& vector,
                       const std::string nm = "dune.gdt.discretefunction")
{
  return DiscreteFunction<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, R>(
      space, std::move(vector.as_imp()), nm);
}


template <class VectorType, class GV, size_t r, size_t rC, class R>
DiscreteFunction<VectorType, GV, r, rC, R> make_discrete_function(const SpaceInterface<GV, r, rC, R>& space,
                                                                  const std::string nm = "dune.gdt.discretefunction")
{
  return DiscreteFunction<VectorType, GV, r, rC, R>(space, nm);
}


template <class GV, size_t r, size_t rC, class R>
DiscreteFunction<typename XT::LA::Container<R>::VectorType, GV, r, rC, R>
make_discrete_function(const SpaceInterface<GV, r, rC, R>& space, const std::string nm = "dune.gdt.discretefunction")
{
  return DiscreteFunction<typename XT::LA::Container<R>::VectorType, GV, r, rC, R>(space, nm);
}


} // namespace GDT
} // namespace Dune

#include "default-datahandle.hh"

#endif // DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
