// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2018)

#ifndef DUNE_GDT_MOMENTMODELS_ENTROPYSOLVER_HH
#define DUNE_GDT_MOMENTMODELS_ENTROPYSOLVER_HH

#include <string>

#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/momentmodels/entropybased_flux.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class SpaceType, class VectorType, class BasisfunctionType>
class LocalEntropySolver : public XT::Grid::ElementFunctor<typename SpaceType::GridViewType>
{
  using GridViewType = typename SpaceType::GridViewType;
  using EntityType = typename GridViewType::template Codim<0>::Entity;
  using IndexSetType = typename GridViewType::IndexSet;
  using EntropyFluxType = EntropyBasedLocalFlux<BasisfunctionType>;
  using RangeFieldType = typename EntropyFluxType::RangeFieldType;
  using LocalVectorType = typename EntropyFluxType::VectorType;
  static const size_t dimDomain = EntropyFluxType::basis_dimDomain;
  static const size_t dimRange = EntropyFluxType::basis_dimRange;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;

public:
  explicit LocalEntropySolver(const SpaceType& space,
                              const VectorType& source_dofs,
                              VectorType& range_dofs,
                              const EntropyFluxType& analytical_flux,
                              std::vector<LocalVectorType>& alphas,
                              const RangeFieldType min_acceptable_density,
                              const XT::Common::Parameter& param,
                              const std::string filename = "")
    : space_(space)
    , source_(space_, source_dofs, "source")
    , range_(space_, range_dofs, "range")
    , analytical_flux_(analytical_flux)
    , alphas_(alphas)
    , min_acceptable_density_(min_acceptable_density)
    , param_(param)
    , filename_(filename.empty() ? "" : filename + "_regularization.txt")
    , index_set_(space_.grid_view().indexSet())
  {}

  virtual XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalEntropySolver(space_,
                                  source_.dofs().vector(),
                                  range_.dofs().vector(),
                                  analytical_flux_,
                                  alphas_,
                                  min_acceptable_density_,
                                  param_,
                                  filename_);
  }

  void apply_local(const EntityType& entity) override final
  {
    const auto local_source = source_.local_discrete_function(entity);
    auto local_range = range_.local_discrete_function(entity);
    XT::Common::FieldVector<RangeFieldType, dimRange> u;
    for (size_t ii = 0; ii < dimRange; ++ii)
      u[ii] = local_source->dofs().get_entry(ii);
    const auto& basis_functions = analytical_flux_.basis_functions();
    thread_local auto vector_indices = source_.space().mapper().global_indices(entity);
    source_.space().mapper().global_indices(entity, vector_indices);
    basis_functions.ensure_min_density(u, min_acceptable_density_);
    const auto alpha_ret = analytical_flux_.get_alpha(u, param_, true);
    alphas_[index_set_.index(entity)] = alpha_ret->first;
    const auto regularization_params = alpha_ret->second;
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_range->dofs().set_entry(ii, regularization_params.first[ii]);
    const auto s = regularization_params.second;
    if (s > 0.) {
      if (!filename_.empty()) {
        static std::mutex outfile_lock;
        outfile_lock.lock();
        std::ofstream outfile(filename_, std::ios_base::app);
        outfile << param_.get("t")[0];
        for (size_t ii = 0; ii < dimDomain; ++ii)
          outfile << " " << entity.geometry().center()[ii];
        outfile << " " << s << " ";
        outfile << XT::Common::to_string(u, 15) << std::endl;
        outfile_lock.unlock();
      }
    }
  } // void apply_local(...)

private:
  const SpaceType& space_;
  const ConstDiscreteFunctionType source_;
  DiscreteFunctionType range_;
  const EntropyFluxType& analytical_flux_;
  std::vector<LocalVectorType>& alphas_;
  const RangeFieldType min_acceptable_density_;
  const XT::Common::Parameter& param_;
  const std::string filename_;
  const typename SpaceType::GridViewType::IndexSet& index_set_;
}; // class LocalEntropySolver<...>


template <class BasisfunctionImp,
          class SpaceImp,
          class MatrixType = typename XT::LA::Container<typename BasisfunctionImp::RangeFieldType>::MatrixType>
class EntropySolver
  : public OperatorInterface<MatrixType, typename SpaceImp::GridViewType, BasisfunctionImp::dimRange, 1>
{
  using BaseType = OperatorInterface<MatrixType, typename SpaceImp::GridViewType, BasisfunctionImp::dimRange, 1>;

public:
  using typename BaseType::VectorType;
  using BasisfunctionType = BasisfunctionImp;
  using SpaceType = SpaceImp;
  using SourceSpaceType = SpaceImp;
  using RangeSpaceType = SpaceImp;
  using EntropyFluxType = EntropyBasedLocalFlux<BasisfunctionType>;
  using RangeFieldType = typename BasisfunctionType::RangeFieldType;
  using LocalVectorType = typename EntropyFluxType::VectorType;

  EntropySolver(const EntropyFluxType& analytical_flux,
                const SpaceType& space,
                std::vector<LocalVectorType>& alphas,
                const RangeFieldType min_acceptable_density,
                const std::string filename = "")
    : analytical_flux_(analytical_flux)
    , space_(space)
    , alphas_(alphas)
    , min_acceptable_density_(min_acceptable_density)
    , filename_(filename)
  {}

  virtual bool linear() const override final
  {
    return false;
  }

  virtual const SourceSpaceType& source_space() const override final
  {
    return space_;
  }

  virtual const RangeSpaceType& range_space() const override final
  {
    return space_;
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    LocalEntropySolver<SpaceType, VectorType, BasisfunctionType> local_entropy_solver(
        space_, source, range, analytical_flux_, alphas_, min_acceptable_density_, param, filename_);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_entropy_solver);
    walker.walk(true);
  } // void apply(...)

private:
  const EntropyFluxType& analytical_flux_;
  const SpaceType& space_;
  std::vector<LocalVectorType>& alphas_;
  const RangeFieldType min_acceptable_density_;
  const std::string filename_;
}; // class EntropySolver<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MOMENTMODELS_ENTROPYSOLVER_HH
