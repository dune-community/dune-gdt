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
  typedef typename SourceType::SpaceType SpaceType;
  typedef typename SpaceType::GridLayerType GridLayerType;
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  using EntropyFluxType = EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, SourceType>;

public:
  explicit LocalMomentRegularizer(const SourceType& source,
                                  RangeType& range,
                                  const AnalyticalFluxType& analytical_flux,
                                  const XT::Common::Parameter& param,
                                  const std::string filename)
    : source_(source)
    , range_(range)
    , analytical_flux_(analytical_flux)
    , param_(param)
    , filename_(filename)
  {
    param_.set("boundary", {0.});
  }

  void apply_local(const EntityType& entity)
  {
    const auto x_in_inside_coords = entity.geometry().local(entity.geometry().center());
    const auto u = source_.local_function(entity)->evaluate(x_in_inside_coords, param_);
    assert(dynamic_cast<const EntropyFluxType*>(&analytical_flux_) != nullptr
           && "analytical_flux_ has to be derived from EntropyBasedLocalFlux");
    const auto s = dynamic_cast<const EntropyFluxType*>(&analytical_flux_)
                       ->derived_local_function(entity)
                       ->get_alpha(x_in_inside_coords, u, param_, true, false)
                       .second;

    thread_local auto vector_indices = source_.space().mapper().globalIndices(entity);
    source_.space().mapper().globalIndices(entity, vector_indices);

    // if regularization was needed, we also need to replace u_n in that cell by its regularized version
    if (s > 0.) {
      if (!filename_.empty()) {
        static std::mutex outfile_lock;
        outfile_lock.lock();
        std::ofstream outfile(filename_, std::ios_base::app);
        outfile << param_.get("t")[0];
        for (size_t ii = 0; ii < AnalyticalFluxType::dimDomain; ++ii)
          outfile << " " << entity.geometry().center()[ii];
        outfile << " " << s << " 0" << std::endl;
        outfile_lock.unlock();
      }
      const auto& basis_functions = dynamic_cast<const EntropyFluxType*>(&analytical_flux_)->basis_functions();
      const auto u_iso_scaled = basis_functions.u_iso() * basis_functions.density(u);
      for (size_t ii = 0; ii < vector_indices.size(); ++ii)
        range_.vector().set_entry(vector_indices[ii],
                                  (1 - s) * source_.vector().get_entry(vector_indices[ii]) + s * u_iso_scaled[ii]);
    } else {
      for (const auto& index : vector_indices)
        range_.vector().set_entry(index, source_.vector().get_entry(index));
    }
  } // void apply_local(...)

private:
  const SourceType& source_;
  RangeType& range_;
  const AnalyticalFluxType& analytical_flux_;
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

  MomentRegularizer(const AnalyticalFluxType& analytical_flux, const std::string filename = "")
    : analytical_flux_(analytical_flux)
    , filename_(filename)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    LocalMomentRegularizer<SourceType, RangeType, AnalyticalFluxType, BasisfunctionType> local_moment_regularizer(
        source, range, analytical_flux_, param, filename_);
    auto walker = XT::Grid::Walker<typename SourceType::SpaceType::GridLayerType>(source.space().grid_layer());
    walker.append(local_moment_regularizer);
    walker.walk(true);
  } // void apply(...)

private:
  const AnalyticalFluxType& analytical_flux_;
  const std::string filename_;
}; // class MomentRegularizer<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_REGULARIZATION_HH
