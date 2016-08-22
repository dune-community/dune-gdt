// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_BASE_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_BASE_HH

#include <dune/xt/common/memory.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class ProblemBase : public ProblemInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
private:
  typedef ProblemInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;

public:
  using typename BaseType::DiffusionFactorType;
  using typename BaseType::DiffusionTensorType;
  using typename BaseType::FunctionType;

  ProblemBase(const DiffusionFactorType& diff_fac, const DiffusionTensorType& diff_ten, const FunctionType& forc,
              const FunctionType& dir, const FunctionType& neum, XT::Common::Configuration grd_cfg,
              XT::Common::Configuration bnd_cfg)
    : diffusion_factor_(diff_fac)
    , diffusion_tensor_(diff_ten)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neum)
    , grid_cfg_(grd_cfg)
    , boundary_info_cfg_(bnd_cfg)
  {
  }

  /**
   * \note Do not manually delete these pointers, they are managed automaticall from here on!
   */
  ProblemBase(const DiffusionFactorType* diff_fac, const DiffusionTensorType* diff_ten, const FunctionType* forc,
              const FunctionType* dir, const FunctionType* neum, XT::Common::Configuration grd_cfg,
              XT::Common::Configuration bnd_cfg)
    : diffusion_factor_(diff_fac)
    , diffusion_tensor_(diff_ten)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neum)
    , grid_cfg_(grd_cfg)
    , boundary_info_cfg_(bnd_cfg)
  {
  }

  virtual const DiffusionFactorType& diffusion_factor() const override
  {
    return diffusion_factor_.access();
  }

  virtual const DiffusionTensorType& diffusion_tensor() const override
  {
    return diffusion_tensor_.access();
  }

  virtual const FunctionType& force() const override
  {
    return force_.access();
  }

  virtual const FunctionType& dirichlet() const override
  {
    return dirichlet_.access();
  }

  virtual const FunctionType& neumann() const override
  {
    return neumann_.access();
  }

  virtual const XT::Common::Configuration& grid_cfg() const override
  {
    return grid_cfg_;
  }

  virtual const XT::Common::Configuration& boundary_info_cfg() const override
  {
    return boundary_info_cfg_;
  }

protected:
  const XT::Common::ConstStorageProvider<DiffusionFactorType> diffusion_factor_;
  const XT::Common::ConstStorageProvider<DiffusionTensorType> diffusion_tensor_;
  const XT::Common::ConstStorageProvider<FunctionType> force_;
  const XT::Common::ConstStorageProvider<FunctionType> dirichlet_;
  const XT::Common::ConstStorageProvider<FunctionType> neumann_;
  const XT::Common::Configuration grid_cfg_;
  const XT::Common::Configuration boundary_info_cfg_;
}; // class ProblemBase


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_BASE_HH
