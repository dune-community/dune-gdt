// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_TEST_CELLMODEL_PRECONDITIONERS_HH
#define DUNE_GDT_TEST_CELLMODEL_PRECONDITIONERS_HH

#include <dune/istl/preconditioner.hh>

namespace Dune {


template <class VectorType>
class IdentityPreconditioner : public Dune::Preconditioner<VectorType, VectorType>
{
public:
  using domain_type = VectorType;
  using range_type = VectorType;
  using field_type = typename VectorType::ScalarType;

  IdentityPreconditioner(const SolverCategory::Category cat)
    : category_(cat)
  {}

  //! Category of the preconditioner (see SolverCategory::Category)
  SolverCategory::Category category() const override final
  {
    return category_;
  }

  void pre(domain_type&, range_type&) override final {}

  void apply(domain_type& v, const range_type& d) override final
  {
    v = d;
  }

  void post(domain_type&) override final {}

private:
  SolverCategory::Category category_;
};

template <class VectorType>
class IterativeSolverPreconditioner : public Dune::Preconditioner<VectorType, VectorType>
{
public:
  using domain_type = VectorType;
  using range_type = VectorType;
  using field_type = typename VectorType::ScalarType;
  using SolverType = Dune::InverseOperator<domain_type, range_type>;

  IterativeSolverPreconditioner(std::shared_ptr<SolverType> solver, const SolverCategory::Category cat)
    : solver_(solver)
    , category_(cat)
  {}

  //! Category of the preconditioner (see SolverCategory::Category)
  SolverCategory::Category category() const override final
  {
    return category_;
  }

  void pre(domain_type&, range_type&) override final {}

  void apply(domain_type& v, const range_type& d) override final
  {
    Dune::InverseOperatorResult res;
    auto d2 = d;
    solver_->apply(v, d2, res);
  }

  void post(domain_type&) override final {}

private:
  std::shared_ptr<SolverType> solver_;
  SolverCategory::Category category_;
};


} // namespace Dune

#endif // DUNE_GDT_TEST_CELLMODEL_PRECONDITIONERS_HH
