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
#include <dune/gdt/test/momentmodels/entropyflux.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class SpaceType, class VectorType, class MomentBasis>
class LocalEntropySolver : public XT::Grid::ElementFunctor<typename SpaceType::GridViewType>
{
  using GridViewType = typename SpaceType::GridViewType;
  using BaseType = XT::Grid::ElementFunctor<GridViewType>;
  using EntityType = typename GridViewType::template Codim<0>::Entity;
  using IndexSetType = typename GridViewType::IndexSet;
  using EntropyFluxType = EntropyBasedFluxFunction<GridViewType, MomentBasis>;
  using RangeFieldType = typename EntropyFluxType::RangeFieldType;
  using LocalVectorType = typename EntropyFluxType::VectorType;
  static constexpr size_t dimFlux = EntropyFluxType::dimFlux;
  static constexpr size_t dimRange = EntropyFluxType::basis_dimRange;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;

public:
  explicit LocalEntropySolver(const SpaceType& space,
                              const VectorType& source_dofs,
                              VectorType& range_dofs,
                              const EntropyFluxType& analytical_flux,
                              const RangeFieldType min_acceptable_density,
                              const XT::Common::Parameter& param,
                              const std::string filename = "")
    : space_(space)
    , source_(space_, source_dofs, "source")
    , local_source_(source_.local_discrete_function())
    , range_(space_, range_dofs, "range")
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(analytical_flux)
    , local_flux_(analytical_flux_.derived_local_function())
    , min_acceptable_density_(min_acceptable_density)
    , param_(param)
    , filename_(filename.empty() ? "" : (filename + "_regularization.txt"))
    , index_set_(space_.grid_view().indexSet())
  {}

  explicit LocalEntropySolver(LocalEntropySolver& other)
    : BaseType(other)
    , space_(other.space_)
    , source_(space_, other.source_.dofs().vector(), "source")
    , local_source_(source_.local_discrete_function())
    , range_(space_, other.range_.dofs().vector(), "range")
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(other.analytical_flux_)
    , local_flux_(analytical_flux_.derived_local_function())
    , min_acceptable_density_(other.min_acceptable_density_)
    , param_(other.param_)
    , filename_(other.filename_)
    , index_set_(space_.grid_view().indexSet())
  {}

  XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalEntropySolver(*this);
  }

  void apply_local(const EntityType& entity) override final
  {
    local_source_->bind(entity);
    local_range_->bind(entity);
    XT::Common::FieldVector<RangeFieldType, dimRange> u;
    for (size_t ii = 0; ii < dimRange; ++ii)
      u[ii] = local_source_->dofs().get_entry(ii);
    const auto& basis_functions = analytical_flux_.basis_functions();
    basis_functions.ensure_min_density(u, min_acceptable_density_);
    local_flux_->bind(entity);
    const auto alpha_ret = local_flux_->get_alpha(u, true);
    const auto regularization_params = alpha_ret->second;
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_range_->dofs().set_entry(ii, regularization_params.first[ii]);
    const auto s = regularization_params.second;
    if (s > 0.) {
      if (!filename_.empty()) {
        static std::mutex outfile_lock;
        outfile_lock.lock();
        std::ofstream outfile(filename_, std::ios_base::app);
        outfile << param_.get("t")[0];
        for (size_t ii = 0; ii < dimFlux; ++ii)
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
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_source_;
  DiscreteFunctionType range_;
  std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType> local_range_;
  const EntropyFluxType& analytical_flux_;
  std::unique_ptr<typename EntropyFluxType::Localfunction> local_flux_;
  const RangeFieldType min_acceptable_density_;
  const XT::Common::Parameter& param_;
  const std::string filename_;
  const typename SpaceType::GridViewType::IndexSet& index_set_;
}; // class LocalEntropySolver<...>


template <class MomentBasisImp,
          class SpaceImp,
          class MatrixType = typename XT::LA::Container<typename MomentBasisImp::RangeFieldType>::MatrixType>
class EntropySolver : public OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>
{
  using BaseType = OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>;

public:
  using typename BaseType::VectorType;
  using MomentBasis = MomentBasisImp;
  using SpaceType = SpaceImp;
  using SourceSpaceType = SpaceImp;
  using RangeSpaceType = SpaceImp;
  using EntropyFluxType = EntropyBasedFluxFunction<typename SpaceType::GridViewType, MomentBasis>;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using LocalVectorType = typename EntropyFluxType::VectorType;

  EntropySolver(const EntropyFluxType& analytical_flux,
                const SpaceType& space,
                const RangeFieldType min_acceptable_density,
                const std::string filename = "")
    : analytical_flux_(analytical_flux)
    , space_(space)
    , min_acceptable_density_(min_acceptable_density)
    , filename_(filename)
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return space_;
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    LocalEntropySolver<SpaceType, VectorType, MomentBasis> local_entropy_solver(
        space_, source, range, analytical_flux_, min_acceptable_density_, param, filename_);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_entropy_solver);
    walker.walk(true);
  } // void apply(...)

private:
  const EntropyFluxType& analytical_flux_;
  const SpaceType& space_;
  const RangeFieldType min_acceptable_density_;
  const std::string filename_;
}; // class EntropySolver<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MOMENTMODELS_ENTROPYSOLVER_HH
