// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2018)

#ifndef DUNE_GDT_OPERATORS_FV_REGULARIZATION_HH
#define DUNE_GDT_OPERATORS_FV_REGULARIZATION_HH

#include <dune/xt/grid/walker/functors.hh>

#include <dune/gdt/local/fluxes/entropybased.hh>
#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {


template <class SourceType, class RangeType, class AnalyticalFluxType, class BasisfunctionType>
class LocalMomentRegularizer : public XT::Grid::Functor::Codim0<typename SourceType::SpaceType::GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static_assert(is_discrete_function<SourceType>::value, "SourceType has to be derived from DiscreteFunction!");
  using SpaceType = typename SourceType::SpaceType;
  using GridLayerType = typename SpaceType::GridLayerType;
  using EntityType = typename GridLayerType::template Codim<0>::Entity;
  using IndexSetType = typename GridLayerType::IndexSet;
  using EntropyFluxType = EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, SourceType>;
  using RangeFieldType = typename EntropyFluxType::RangeFieldType;

public:
  explicit LocalMomentRegularizer(const SourceType& source,
                                  RangeType& range,
                                  const AnalyticalFluxType& analytical_flux,
                                  const RangeFieldType min_acceptable_density,
                                  const XT::Common::Parameter& param,
                                  const std::string filename)
    : source_(source)
    , range_(range)
    , analytical_flux_(analytical_flux)
    , min_acceptable_density_(min_acceptable_density)
    , param_(param)
    , filename_(filename + "_regularization.txt")
  {
    param_.set("boundary", {0.});
  }

  void apply_local(const EntityType& entity)
  {
    const auto x_in_inside_coords = entity.geometry().local(entity.geometry().center());
    auto u = source_.local_function(entity)->evaluate(x_in_inside_coords, param_);
    assert(dynamic_cast<const EntropyFluxType*>(&analytical_flux_) != nullptr
           && "analytical_flux_ has to be derived from EntropyBasedLocalFlux");
    const auto* entropy_flux = dynamic_cast<const EntropyFluxType*>(&analytical_flux_);
    const auto& basis_functions = entropy_flux->basis_functions();
    const auto density_min = basis_functions.density_min(u);
    thread_local auto vector_indices = source_.space().mapper().globalIndices(entity);
    source_.space().mapper().globalIndices(entity, vector_indices);
    if (density_min < min_acceptable_density_ / basis_functions.density_factor()) {
      std::cerr << "Added small vaccuum density to " << XT::Common::to_string(u, 15) << " on entity "
                << XT::Common::to_string(entity.geometry().center(), 15) << std::endl;
      u += basis_functions.u_iso() * min_acceptable_density_;
    }
    const auto s =
        entropy_flux->derived_local_function(entity)->get_alpha(x_in_inside_coords, u, param_, true, false).second;
    // if regularization was needed, we also need to replace u_n in that cell by its regularized version
    if (s > 0.) {
      if (!filename_.empty()) {
        static std::mutex outfile_lock;
        outfile_lock.lock();
        std::ofstream outfile(filename_, std::ios_base::app);
        outfile << param_.get("t")[0];
        for (size_t ii = 0; ii < AnalyticalFluxType::dimDomain; ++ii)
          outfile << " " << entity.geometry().center()[ii];
        outfile << " " << s << " ";
        outfile << XT::Common::to_string(u, 15) << std::endl;
        outfile_lock.unlock();
      }
      const auto u_iso_scaled = basis_functions.u_iso() * basis_functions.density(u);
      for (size_t ii = 0; ii < vector_indices.size(); ++ii)
        range_.vector().set_entry(vector_indices[ii], (1 - s) * u[ii] + s * u_iso_scaled[ii]);
    } else {
      for (size_t ii = 0; ii < vector_indices.size(); ++ii)
        range_.vector().set_entry(vector_indices[ii], u[ii]);
    }
  } // void apply_local(...)

private:
  const SourceType& source_;
  RangeType& range_;
  const AnalyticalFluxType& analytical_flux_;
  const RangeFieldType min_acceptable_density_;
  XT::Common::Parameter param_;
  const std::string filename_;
}; // class LocalMomentRegularizer<...>


template <class AnalyticalFluxType, class BasisfunctionType, class Traits>
class MomentRegularizer;


namespace internal {


template <class AnalyticalFluxImp, class BasisfunctionImp>
struct MomentRegularizerTraits
{
  using AnalyticalFluxType = AnalyticalFluxImp;
  using BasisfunctionType = BasisfunctionImp;
  using FieldType = typename AnalyticalFluxType::DomainFieldType;
  using JacobianType = NoJacobian;
  using RangeFieldType = typename BasisfunctionType::RangeFieldType;
  using derived_type = MomentRegularizer<AnalyticalFluxType, BasisfunctionType, MomentRegularizerTraits>;
};


} // namespace internal


template <class AnalyticalFluxImp,
          class BasisfunctionImp,
          class Traits = internal::MomentRegularizerTraits<AnalyticalFluxImp, BasisfunctionImp>>
class MomentRegularizer : public OperatorInterface<Traits>
{
public:
  using AnalyticalFluxType = typename Traits::AnalyticalFluxType;
  using BasisfunctionType = typename Traits::BasisfunctionType;
  using RangeFieldType = typename Traits::RangeFieldType;

  MomentRegularizer(const AnalyticalFluxType& analytical_flux,
                    const RangeFieldType min_acceptable_density,
                    const std::string filename = "")
    : analytical_flux_(analytical_flux)
    , min_acceptable_density_(min_acceptable_density)
    , filename_(filename)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    LocalMomentRegularizer<SourceType, RangeType, AnalyticalFluxType, BasisfunctionType> local_moment_regularizer(
        source, range, analytical_flux_, min_acceptable_density_, param, filename_);
    auto walker = XT::Grid::Walker<typename SourceType::SpaceType::GridLayerType>(source.space().grid_layer());
    walker.append(local_moment_regularizer);
    walker.walk(true);
  } // void apply(...)

private:
  const AnalyticalFluxType& analytical_flux_;
  const RangeFieldType min_acceptable_density_;
  const std::string filename_;
}; // class MomentRegularizer<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_REGULARIZATION_HH
