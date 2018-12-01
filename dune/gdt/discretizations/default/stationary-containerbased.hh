// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_DISCRETIZATIONS_DEFAULT_STATIONARY_CONTAINERBASED_HH
#define DUNE_GDT_DISCRETIZATIONS_DEFAULT_STATIONARY_CONTAINERBASED_HH

#include <dune/xt/common/exceptions.hh>

#include "../interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template <class ProblemType,
          class AnsatzSpaceType,
          class MatrixType,
          class VectorType,
          class TestSpaceType = AnsatzSpaceType>
class StationaryContainerBasedDefaultDiscretization;


namespace internal {


template <class ProblemImp, class AnsatzSpaceImp, class MatrixImp, class VectorImp, class TestSpaceImp>
class StationaryContainerBasedDefaultDiscretizationTraits
{
  // no checks of the arguments needed, those are done in the interfaces
public:
  typedef StationaryContainerBasedDefaultDiscretization<ProblemImp, AnsatzSpaceImp, MatrixImp, VectorImp, TestSpaceImp>
      derived_type;
  typedef ProblemImp ProblemType;
  typedef AnsatzSpaceImp AnsatzSpaceType;
  typedef TestSpaceImp TestSpaceType;
  typedef MatrixImp MatrixType;
  typedef VectorImp VectorType;
}; // class StationaryContainerBasedDefaultDiscretizationTraits


} // namespace internal


template <class ProblemImp, class AnsatzSpaceImp, class MatrixImp, class VectorImp, class TestSpaceImp>
class StationaryContainerBasedDefaultDiscretization
  : public ContainerBasedStationaryDiscretizationInterface<
        internal::StationaryContainerBasedDefaultDiscretizationTraits<ProblemImp,
                                                                      AnsatzSpaceImp,
                                                                      MatrixImp,
                                                                      VectorImp,
                                                                      TestSpaceImp>>
{
  typedef ContainerBasedStationaryDiscretizationInterface<
      internal::StationaryContainerBasedDefaultDiscretizationTraits<ProblemImp,
                                                                    AnsatzSpaceImp,
                                                                    MatrixImp,
                                                                    VectorImp,
                                                                    TestSpaceImp>>
      BaseType;
  typedef StationaryContainerBasedDefaultDiscretization<ProblemImp, AnsatzSpaceImp, MatrixImp, VectorImp, TestSpaceImp>
      ThisType;

public:
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::MatrixType;
  using typename BaseType::ProblemType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::VectorType;

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm,
                                                AnsatzSpaceType ansatz_sp,
                                                TestSpaceType test_sp,
                                                MatrixType system_mtrx,
                                                VectorType rhs_vec,
                                                VectorType dirichlet)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(test_sp)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_(dirichlet)
    , has_dirichlet_shift_(true)
  {}

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm,
                                                AnsatzSpaceType ansatz_sp,
                                                MatrixType system_mtrx,
                                                VectorType rhs_vec,
                                                VectorType dirichlet)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(ansatz_space_)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_(dirichlet)
    , has_dirichlet_shift_(true)
  {}

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm,
                                                AnsatzSpaceType ansatz_sp,
                                                TestSpaceType test_sp,
                                                MatrixType system_mtrx,
                                                VectorType rhs_vec)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(test_sp)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_()
    , has_dirichlet_shift_(false)
  {}

  StationaryContainerBasedDefaultDiscretization(const ProblemType& prblm,
                                                AnsatzSpaceType ansatz_sp,
                                                MatrixType system_mtrx,
                                                VectorType rhs_vec)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(ansatz_space_)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_()
    , has_dirichlet_shift_(false)
  {}

  StationaryContainerBasedDefaultDiscretization(ThisType&& /*source*/) = default;

  /// \name Required by StationaryDiscretizationInterface.
  /// \{

  const ProblemType& problem() const
  {
    return problem_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return ansatz_space_;
  }

  const TestSpaceType& test_space() const
  {
    return test_space_;
  }

  /// \}
  /// \name Required by ContainerBasedStationaryDiscretizationInterface.
  /// \{

  const MatrixType& system_matrix() const
  {
    return system_matrix_;
  }

  const VectorType& rhs_vector() const
  {
    return rhs_vector_;
  }

  bool has_dirichlet_shift() const
  {
    return has_dirichlet_shift_;
  }

  const VectorType& dirichlet_shift() const
  {
    if (!has_dirichlet_shift_)
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Do not call dirichlet_shift() if has_dirichlet_shift() is false!");
    return dirichlet_shift_;
  }

  /// \}

private:
  const ProblemType& problem_;
  const AnsatzSpaceType ansatz_space_;
  const TestSpaceType test_space_;
  const MatrixType system_matrix_;
  const VectorType rhs_vector_;
  const VectorType dirichlet_shift_;
  const bool has_dirichlet_shift_;
}; // class StationaryContainerBasedDefaultDiscretization


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_DEFAULT_STATIONARY_CONTAINERBASED_HH
