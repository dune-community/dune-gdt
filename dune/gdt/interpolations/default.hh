// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_INTERPOLATIONS_DEFAULT_HH
#define DUNE_GDT_INTERPOLATIONS_DEFAULT_HH

#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/walker.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


// ### Variants for GridFunctionInterface
template <class GV, size_t r, size_t rC, class R, class V, class IGV>
class DefaultInterpolationElementFunctor : public XT::Grid::ElementFunctor<IGV>
{
  using BaseType = typename XT::Grid::ElementFunctor<IGV>;

public:
  using typename BaseType::E;
  using SourceType = XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>;
  using LocalSourceType = typename SourceType::LocalFunctionType;
  using TargetType = DiscreteFunction<V, GV, r, rC, R>;
  using LocalDofVectorType = typename TargetType::DofVectorType::LocalDofVectorType;
  using TargetBasisType = typename TargetType::SpaceType::GlobalBasisType::LocalizedType;

  DefaultInterpolationElementFunctor(const SourceType& source, TargetType& target)
    : source_(source)
    , target_(target)
    , local_dof_vector_(target.dofs().localize())
    , local_source_(source.local_function())
    , target_basis_(target.space().basis().localize())
  {
    DUNE_THROW_IF(target_.space().type() == SpaceType::raviart_thomas,
                  Exceptions::interpolation_error,
                  "Use the correct one from interpolations/raviart-thomas.hh instead!");
  }

  DefaultInterpolationElementFunctor(const DefaultInterpolationElementFunctor& other)
    : source_(other.source_)
    , target_(other.target_)
    , local_dof_vector_(target_.dofs().localize())
    , local_source_(source_.local_function())
    , target_basis_(target_.space().basis().localize())
  {}

  virtual XT::Grid::ElementFunctor<IGV>* copy() override final
  {
    return new DefaultInterpolationElementFunctor(*this);
  }

  virtual void apply_local(const E& element) override final
  {
    local_source_->bind(element);
    local_dof_vector_.bind(element);
    target_basis_->bind(element);
    target_basis_->interpolate(
        [&](const auto& xx) { return local_source_->evaluate(xx); }, local_source_->order(), local_dof_vector_);
  }

private:
  const SourceType& source_;
  TargetType& target_;
  LocalDofVectorType local_dof_vector_;
  std::unique_ptr<LocalSourceType> local_source_;
  std::unique_ptr<TargetBasisType> target_basis_;
};

/**
 * \brief Interpolates a localizable function within a given space [most general variant].
 *
 * Simply uses interpolate() of the target spaces global basis().
 *
 *
 * \note This does not clear target.dofs().vector(). Thus, if interpolation_grid_view only covers a part of the domain
 *       of target.space().grid_view(), other contributions in target remain (which is on purpose).
 *
 * \note This might not be optimal for all spaces. For instance, the polynomial order of source is not taken into
 *       account for local L^2-projection based interpolation. This is a limitation in dune-localfunctions and we need
 *       to completely replace the interpolation of the respective local finite element.
 *
 * \note This might not be correct for all spaces, in particular if source is not continuous.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGVT>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGVT::Grid::template Codim<0>::Entity>::value,
                 void>
default_interpolation(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
                      DiscreteFunction<V, GV, r, rC, R>& target,
                      const GridView<IGVT>& interpolation_grid_view)
{
  DefaultInterpolationElementFunctor<GV, r, rC, R, V, GridView<IGVT>> functor(source, target);
  auto walker = XT::Grid::Walker<GridView<IGVT>>(interpolation_grid_view);
  walker.append(functor);
  walker.walk(true);
} // ... default_interpolation(...)


/**
 * \brief Interpolates a localizable function within a given space [uses target.space().grid_view() as
 *        interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void default_interpolation(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
                           DiscreteFunction<V, GV, r, rC, R>& target)
{
  default_interpolation(source, target, target.space().grid_view());
}


/**
 * \brief Interpolates a localizable function within a given space [creates a suitable target function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGVT>
std::enable_if_t<
    XT::LA::is_vector<VectorType>::value
        && std::is_same<XT::Grid::extract_entity_t<GV>, typename IGVT::Grid::template Codim<0>::Entity>::value,
    DiscreteFunction<VectorType, GV, r, rC, R>>
default_interpolation(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
                      const SpaceInterface<GV, r, rC, R>& target_space,
                      const GridView<IGVT>& interpolation_grid_view)
{
  auto target_function = make_discrete_function<VectorType>(target_space);
  default_interpolation(source, target_function, interpolation_grid_view);
  return target_function;
}


/**
 * \brief Interpolates a localizable function within a given space [creates a suitable target function, uses
 *        target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
default_interpolation(const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source,
                      const SpaceInterface<GV, r, rC, R>& target_space)
{
  auto target_function = make_discrete_function<VectorType>(target_space);
  default_interpolation(source, target_function);
  return target_function;
}


// ### Variants for FunctionInterface


/**
 * \brief Interpolates a function within a given space [most general variant].
 *
 * Simply calls as_grid_function<>() and redirects to the appropriate default_interpolation() function.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGVT>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGVT::Grid::template Codim<0>::Entity>::value,
                 void>
default_interpolation(const XT::Functions::FunctionInterface<GridView<IGVT>::dimension, r, rC, R>& source,
                      DiscreteFunction<V, GV, r, rC, R>& target,
                      const GridView<IGVT>& interpolation_grid_view)
{
  default_interpolation(source.as_grid_function(interpolation_grid_view), target, interpolation_grid_view);
}


/**
 * \brief Interpolates a function within a given space [uses target.space().grid_view() as
 *        interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void default_interpolation(const XT::Functions::FunctionInterface<GV::dimension, r, rC, R>& source,
                           DiscreteFunction<V, GV, r, rC, R>& target)
{
  default_interpolation(source, target, target.space().grid_view());
}


/**
 * \brief Interpolates a function within a given space [creates a suitable target function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGVT>
std::enable_if_t<
    XT::LA::is_vector<VectorType>::value
        && std::is_same<XT::Grid::extract_entity_t<GV>, typename IGVT::Grid::template Codim<0>::Entity>::value,
    DiscreteFunction<VectorType, GV, r, rC, R>>
default_interpolation(const XT::Functions::FunctionInterface<GridView<IGVT>::dimension, r, rC, R>& source,
                      const SpaceInterface<GV, r, rC, R>& target_space,
                      const GridView<IGVT>& interpolation_grid_view)
{
  return default_interpolation<VectorType>(
      source.as_grid_function(interpolation_grid_view), target_space, interpolation_grid_view);
}


/**
 * \brief Interpolates a function within a given space [creates a suitable target function, uses
 *        target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
default_interpolation(const XT::Functions::FunctionInterface<GV::dimension, r, rC, R>& source,
                      const SpaceInterface<GV, r, rC, R>& target_space)
{
  return default_interpolation<VectorType>(source, target_space, target_space.grid_view());
}


// ### Variants for GenericFunction


/**
 * \brief Interpolates a function given as a lambda expression within a given space [most general variant].
 *
 * Simply creates a XT::Functions::GenericFunction and redirects to the appropriate default_interpolation() function.
 */
template <class GV, size_t r, size_t rC, class R, class V, class IGVT>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<GV>, typename IGVT::Grid::template Codim<0>::Entity>::value,
                 void>
default_interpolation(
    const int source_order,
    const std::function<typename XT::Functions::GenericFunction<GridView<IGVT>::dimension, r, rC, R>::RangeReturnType(
        const typename XT::Functions::GenericFunction<GridView<IGVT>::dimension, r, rC, R>::DomainType&,
        const XT::Common::Parameter&)> source_evaluate_lambda,
    DiscreteFunction<V, GV, r, rC, R>& target,
    const GridView<IGVT>& interpolation_grid_view)
{
  default_interpolation(
      XT::Functions::GenericFunction<GridView<IGVT>::dimension, r, rC, R>(source_order, source_evaluate_lambda),
      target,
      interpolation_grid_view);
}


/**
 * \brief Interpolates a function given as a lambda expression within a given space [uses
 *        target.space().grid_view() as interpolation_grid_view].
 **/
template <class GV, size_t r, size_t rC, class R, class V>
void default_interpolation(
    const int source_order,
    const std::function<typename XT::Functions::GenericFunction<GV::dimension, r, rC, R>::RangeReturnType(
        const typename XT::Functions::GenericFunction<GV::dimension, r, rC, R>::DomainType&,
        const XT::Common::Parameter&)> source_evaluate_lambda,
    DiscreteFunction<V, GV, r, rC, R>& target)
{
  default_interpolation(XT::Functions::GenericFunction<GV::dimension, r, rC, R>(source_order, source_evaluate_lambda),
                        target);
}


/**
 * \brief Interpolates a function given as a lambda expression within a given space [creates a suitable target
 *        function].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R, class IGVT>
std::enable_if_t<
    XT::LA::is_vector<VectorType>::value
        && std::is_same<XT::Grid::extract_entity_t<GV>, typename IGVT::Grid::template Codim<0>::Entity>::value,
    DiscreteFunction<VectorType, GV, r, rC, R>>
default_interpolation(
    const int source_order,
    const std::function<typename XT::Functions::GenericFunction<GridView<IGVT>::dimension, r, rC, R>::RangeReturnType(
        const typename XT::Functions::GenericFunction<GridView<IGVT>::dimension, r, rC, R>::DomainType&,
        const XT::Common::Parameter&)> source_evaluate_lambda,
    const SpaceInterface<GV, r, rC, R>& target_space,
    const GridView<IGVT>& interpolation_grid_view)
{
  return default_interpolation<VectorType>(
      XT::Functions::GenericFunction<GridView<IGVT>::dimension, r, rC, R>(source_order, source_evaluate_lambda),
      target_space,
      interpolation_grid_view);
}


/**
 * \brief Interpolates a function given as a lambda expression within a given space [creates a suitable target
 *        function, uses target_space.grid_view() as interpolation_grid_view].
 **/
template <class VectorType, class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::LA::is_vector<VectorType>::value, DiscreteFunction<VectorType, GV, r, rC, R>>
default_interpolation(
    const int source_order,
    const std::function<typename XT::Functions::GenericFunction<GV::dimension, r, rC, R>::RangeReturnType(
        const typename XT::Functions::GenericFunction<GV::dimension, r, rC, R>::DomainType&,
        const XT::Common::Parameter&)> source_evaluate_lambda,
    const SpaceInterface<GV, r, rC, R>& target_space)
{
  return default_interpolation<VectorType>(
      XT::Functions::GenericFunction<GV::dimension, r, rC, R>(source_order, source_evaluate_lambda), target_space);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_INTERPOLATIONS_DEFAULT_HH
