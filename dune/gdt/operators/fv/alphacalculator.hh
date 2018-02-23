// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_OPERATORS_FV_MOMENTREGULARIZER_HH
#define DUNE_GDT_OPERATORS_FV_MOMENTREGULARIZER_HH

#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/walker/functors.hh>

#include <dune/xt/la/algorithms/triangular_solves.hh>
#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/eigen-solver.hh>

#include "slopelimiters.hh"

namespace Dune {
namespace GDT {


template <class SourceType, class AnalyticalFluxType>
class MomentRegularizer : public XT::Grid::Functor::Codim0<typename SourceType::SpaceType::GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  typedef typename SourceType::SpaceType SpaceType;
  typedef typename SpaceType::GridLayerType GridLayerType;
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;

public:
  explicit MomentRegularizer(SourceType& source,
                             const AnalyticalFluxType& analytical_flux,
                             const XT::Common::Parameter& param)
    : source_(source)
    , analytical_flux_(analytical_flux)
    , grid_layer_(source_.space().grid_layer())
    , param_(param)
  {
    param_.set("boundary", {0.});
  }

  void apply_local(const EntityType& entity)
  {
    const auto x_in_inside_coords = entity.geometry().local(entity.geometry().center());
    const auto u = source_.local_function(entity)->evaluate(x_in_inside_coords, param_);
    analytical_flux_.local_function(entity)->calculate_alpha(x_in_inside_coords, u, param_);

    // if regularization was needed, we also need to replace u_n in that cell by its regularized version
    const auto& entity_index = grid_layer_.indexSet().index(entity);
    const auto& s = analytical_flux_.regularization_parameters()[entity_index];
    if (s > 0.) {
      std::cout << "regularization in regularizer: time: " << param_.get("t")[0]
                << ", entitycenter: " << XT::Common::to_string(entity.geometry().center())
                << ", s: " << XT::Common::to_string(s, 15) << std::endl;
      const auto& vector_indices = source_.space().mapper().globalIndices(entity);
      const auto u_iso = analytical_flux_.get_isotropic_moment(u);
      auto& source_vector = source_.vector();
      for (size_t ii = 0; ii < vector_indices.size(); ++ii) {
        const auto& vec_index = vector_indices[ii];
        source_vector.set_entry(vec_index, (1 - s) * source_vector.get_entry(vec_index) + s * u_iso[ii]);
      }
    }
  } // void apply_local(...)

private:
  SourceType& source_;
  const AnalyticalFluxType& analytical_flux_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
}; // class MomentRegularizer<...>

template <class AnalyticalFluxType, class DiscreteFunctionType, size_t dimDomain, size_t dimRange>
class EntropyProblemSolver : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>
{
  typedef typename DiscreteFunctionType::SpaceType::GridLayerType GridLayerType;
  typedef typename DiscreteFunctionType::EntityType EntityType;
  typedef typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate DomainType;
  typedef typename DiscreteFunctionType::RangeType RangeType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename Dune::QuadratureRule<RangeFieldType, dimDomain> QuadratureType;
  typedef typename std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>> ReconstructedValuesType;

public:
  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit EntropyProblemSolver(
      const AnalyticalFluxType& analytical_flux,
      const std::vector<RangeType>& source_values,
      const IndexSetType& index_set,
      ReconstructedValuesType& reconstructed_values,
      const XT::Common::Parameter param,
      const std::vector<RangeFieldType>& r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1})
    : analytical_flux_(analytical_flux)
    , source_values_(source_values)
    , index_set_(index_set)
    , reconstructed_values_(reconstructed_values)
    , param_(param)
    , r_sequence_(r_sequence)
  {
    param_.set("boundary", {0.});
  }

  void apply_local(const EntityType& entity)
  {
    const auto entity_index = index_set_.index(entity);
    auto& local_reconstructed_values = reconstructed_values_[entity_index];
    const RangeType& u_bar = source_values_[entity_index];
    auto local_func = analytical_flux_.local_function(entity);

    auto limited_values = local_reconstructed_values;
    bool finished = false;
    const auto x_in_inside_coords = entity.geometry().local(entity.geometry().center());
    const auto r_max = r_sequence_.back();
    size_t solves = 0;
    while (!finished) {
      finished = true;
      for (auto& pair : limited_values) {
        auto& u = pair.second;
        RangeFieldType r_used = 0.;
        for (const auto& r : r_sequence_) {
          r_used = r;
          RangeType u_limited = convex_combination(u, u_bar, r);
          try {
            local_func->calculate_alpha(x_in_inside_coords, u_limited, param_);
            ++solves;
            break;
          } catch (const Dune::MathError& e) {
            if (r < r_max)
              continue;
            else
              DUNE_THROW(Dune::MathError, "This was the original error:\n" + XT::Common::to_string(e.what()));
          } // try ... catch ...
        } // r
        if (r_used > 0.) {
          std::cout << "limited in solver with r = " << XT::Common::to_string(r_used)
                    << ", entity: " << XT::Common::to_string(entity.geometry().center()) << std::endl;
          limit_realizability(limited_values, local_reconstructed_values, u_bar, r_used);
          finished = false;
          break;
        }
      } // reconstructed values
    } // while (!finished)
    if (solves > 2 * dimDomain)
      std::cout << "Took: " << solves << std::endl;
    local_reconstructed_values = limited_values;
  } // void apply_local(...)

private:
  void limit_realizability(std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& limited_values,
                           const std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& old_values,
                           const RangeType& u_bar,
                           const RangeFieldType theta)
  {
    for (auto& pair : limited_values) {
      auto& u = pair.second;
      const auto& u_old = old_values.at(pair.first);
      u = convex_combination(u_old, u_bar, theta);
    }
  }

  RangeType convex_combination(const RangeType& u, const RangeType& u_bar, const RangeFieldType& theta)
  {
    if (XT::Common::FloatCmp::eq(theta, 0.))
      return u;
    else {
      RangeType u_scaled = u;
      u_scaled *= 1 - theta;
      RangeType u_bar_scaled = u_bar;
      u_bar_scaled *= theta;
      return u_scaled + u_bar_scaled;
    }
  }

  const AnalyticalFluxType& analytical_flux_;
  const std::vector<RangeType> source_values_;
  const IndexSetType& index_set_;
  ReconstructedValuesType& reconstructed_values_;
  XT::Common::Parameter param_;
  const std::vector<RangeFieldType> r_sequence_;
}; // class EntropyProblemSolver


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_MOMENTREGULARIZER_HH
