// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_DISCRETIZATIONS_INTERFACES_HH
#define DUNE_GDT_DISCRETIZATIONS_INTERFACES_HH

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


template <class Traits>
class StationaryDiscretizationInterface : public Stuff::CRTPInterface<StationaryDiscretizationInterface<Traits>, Traits>
{
  typedef Stuff::CRTPInterface<StationaryDiscretizationInterface<Traits>, Traits> BaseType;

public:
  using typename BaseType::derived_type;
  typedef typename Traits::ProblemType ProblemType;
  typedef typename Traits::AnsatzSpaceType AnsatzSpaceType;
  typedef typename Traits::TestSpaceType TestSpaceType;
  typedef typename Traits::VectorType VectorType;

private:
  static_assert(is_space<AnsatzSpaceType>::value, "AnsatzSpaceType has to be derived from SpaceInterface!");
  static_assert(is_space<TestSpaceType>::value, "TestSpaceType has to be derived from SpaceInterface!");
  static_assert(Stuff::LA::is_vector<VectorType>::value,
                "VectorType has to be derived from Stuff::LA::VectorInterface!");

public:
  /// \name Have to be implemented by any derived class.
  /// \{

  const ProblemType& problem() const
  {
    CHECK_CRTP(this->as_imp().problem());
    return this->as_imp().problem();
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    CHECK_CRTP(this->as_imp().ansatz_space());
    return this->as_imp().ansatz_space();
  }

  const TestSpaceType& test_space() const
  {
    CHECK_CRTP(this->as_imp().test_space());
    return this->as_imp().test_space();
  }

  std::vector<std::string> solver_types() const
  {
    CHECK_CRTP(this->as_imp().solver_types());
    auto types = this->as_imp().solver_types();
    if (types.empty())
      DUNE_THROW(Stuff::Exceptions::internal_error,
                 "Reported solver_types() of the derived class (see below) must not be empty!\n\n  "
                     << Stuff::Common::Typename<derived_type>::value());
    return types;
  } // ... solver_types(...)

  Stuff::Common::Configuration solver_options(const std::string type = "") const
  {
    CHECK_CRTP(this->as_imp().solver_options(type));
    auto opts = this->as_imp().solver_options(type);
    if (opts.empty())
      DUNE_THROW(Stuff::Exceptions::internal_error,
                 "Reported solver_options() of the derived class (see below) for type '"
                     << type
                     << "'must not be empty!\n\n  "
                     << Stuff::Common::Typename<derived_type>::value());
    return opts;
  } // ... solver_options(...)

  void solve(VectorType& solution, const Stuff::Common::Configuration& options) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().solve(solution, options));
  }

  /// \}
  /// \name Provided by the interface for convenience.
  /// \{

  VectorType create_vector() const
  {
    return VectorType(ansatz_space().mapper().size());
  }

  void solve(VectorType& solution, const std::string& type) const
  {
    solve(solution, solver_options(type));
  }

  void solve(VectorType& solution) const
  {
    solve(solution, solver_options(solver_types().at(0)));
  }

  VectorType solve(const Stuff::Common::Configuration& options) const
  {
    VectorType solution = create_vector();
    solve(solution, options);
    return solution;
  }

  VectorType solve(const std::string& type) const
  {
    VectorType solution = create_vector();
    solve(solution, type);
    return solution;
  }

  VectorType solve() const
  {
    VectorType solution = create_vector();
    solve(solution);
    return solution;
  }

  void visualize(const VectorType& vector, const std::string filename, const std::string name) const
  {
    make_const_discrete_function(this->ansatz_space(), vector, name).visualize(filename);
  }

  /// \}
}; // class StationaryDiscretizationInterface


template <class Traits>
class ContainerBasedStationaryDiscretizationInterface : public StationaryDiscretizationInterface<Traits>
{
  typedef StationaryDiscretizationInterface<Traits> BaseType;

public:
  typedef typename Traits::MatrixType MatrixType;
  using typename BaseType::VectorType;

private:
  typedef typename Stuff::LA::Solver<MatrixType> LinearSolverType;

public:
  /// \name Have to be implemented by any derived class.
  /// \{

  const MatrixType& system_matrix() const
  {
    CHECK_CRTP(this->as_imp().system_matrix());
    return this->as_imp().system_matrix();
  }

  const VectorType& rhs_vector() const
  {
    CHECK_CRTP(this->as_imp().rhs_vector());
    return this->as_imp().rhs_vector();
  }

  bool has_dirichlet_shift() const
  {
    CHECK_CRTP(this->as_imp().has_dirichlet_shift());
    return this->as_imp().has_dirichlet_shift();
  }

  /// \}

  /**
   * \brief Returns the Dirichlet shift.
   * \note  This method has to be implemented, if has_dirichlet_shift() returns false!
   */
  const VectorType& dirichlet_shift() const
  {
    if (!has_dirichlet_shift())
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Do not call dirichlet_shift() if has_dirichlet_shift() is false!");
    CHECK_CRTP(this->as_imp().dirichlet_shift());
    return this->as_imp().dirichlet_shift();
  }

  /// \name Provided by the interface for convenience.
  /// \{

  std::vector<std::string> solver_types() const
  {
    return LinearSolverType::types();
  }

  Stuff::Common::Configuration solver_options(const std::string type = "") const
  {
    return LinearSolverType::options(type);
  }

  using BaseType::solve;

  void solve(VectorType& solution, const Stuff::Common::Configuration& options) const
  {
    LinearSolverType(system_matrix()).apply(rhs_vector(), solution, options);
    if (has_dirichlet_shift())
      solution += dirichlet_shift();
  }

  /// \}
}; // class ContainerBasedStationaryDiscretizationInterface


template <class Traits>
class FVDiscretizationInterface : public Stuff::CRTPInterface<FVDiscretizationInterface<Traits>, Traits>
{
  typedef Stuff::CRTPInterface<FVDiscretizationInterface<Traits>, Traits> BaseType;

public:
  using typename BaseType::derived_type;
  typedef typename Traits::ProblemType ProblemType;
  typedef typename Traits::SpaceType SpaceType;
  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::DiscreteSolutionType DiscreteSolutionType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;

private:
  static_assert(is_space<SpaceType>::value, "SpaceType has to be derived from SpaceInterface!");

public:
  /// \name Have to be implemented by any derived class.
  /// \{

  const ProblemType& problem() const
  {
    CHECK_CRTP(this->as_imp().problem());
    return this->as_imp().problem();
  }

  const SpaceType& space() const
  {
    CHECK_CRTP(this->as_imp().space());
    return this->as_imp().space();
  }

  void solve(DiscreteSolutionType& solution) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().solve(solution));
  }

  /// \}
  /// \name Provided by the interface for convenience.
  /// \{

  DiscreteSolutionType solve() const
  {
    DiscreteSolutionType solution;
    solve(solution);
    return solution;
  }

  /// \}
}; // class FVDiscretizationInterface


namespace internal {


template <class D>
struct is_stationary_discretization_helper
{
  DSC_has_typedef_initialize_once(Traits)

      static const bool is_candidate = DSC_has_typedef(Traits)<D>::value;
};


} // namespace internal


template <class D, bool candidate = internal::is_stationary_discretization_helper<D>::is_candidate>
struct is_stationary_discretization : public std::is_base_of<StationaryDiscretizationInterface<typename D::Traits>, D>
{
};


template <class D>
struct is_stationary_discretization<D, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_INTERFACES_HH
