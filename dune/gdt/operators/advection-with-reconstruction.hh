// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_WITH_RECONSTRUCTION_HH
#define DUNE_GDT_OPERATORS_ADVECTION_WITH_RECONSTRUCTION_HH

#include <dune/gdt/operators/interfaces.hh>

#include "reconstruction/linear.hh"

namespace Dune {
namespace GDT {


template <class AdvectionOperatorImp, class ReconstructionOperatorImp>
class AdvectionWithReconstructionOperator
  : public OperatorInterface<typename AdvectionOperatorImp::AGV,
                             AdvectionOperatorImp::s_r,
                             AdvectionOperatorImp::s_rC,
                             AdvectionOperatorImp::r_r,
                             AdvectionOperatorImp::r_rC,
                             typename AdvectionOperatorImp::F,
                             typename AdvectionOperatorImp::M,
                             typename AdvectionOperatorImp::SGV,
                             typename AdvectionOperatorImp::RGV>
{
public:
  using ThisType = AdvectionWithReconstructionOperator;
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::AGV,
                                     AdvectionOperatorImp::s_r,
                                     AdvectionOperatorImp::s_rC,
                                     AdvectionOperatorImp::r_r,
                                     AdvectionOperatorImp::r_rC,
                                     typename AdvectionOperatorImp::F,
                                     typename AdvectionOperatorImp::M,
                                     typename AdvectionOperatorImp::SGV,
                                     typename AdvectionOperatorImp::RGV>;

  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AdvectionOperatorType = AdvectionOperatorImp;
  using ReconstructionOperatorType = ReconstructionOperatorImp;

  AdvectionWithReconstructionOperator(const AdvectionOperatorType& advection_operator,
                                      const ReconstructionOperatorType& reconstruction_operator)
    : BaseType(advection_operator.parameter_type() + reconstruction_operator.parameter_type())
    , advection_operator_(advection_operator)
    , reconstruction_operator_(reconstruction_operator)
    , reconstruction_(reconstruction_operator.range_space().mapper().size())
  {}

  // pull in methods from various base classes
  using BaseType::apply;

  /// \name Required by ForwardOperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override final
  {
    return advection_operator_.range_space();
  }

  bool linear() const override final
  {
    return false;
  }

  // avoid non-optimal default implementation in OperatorInterface
  void apply(SourceFunctionType source_function,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_range(range_vector);
    reconstruction_operator_.apply(source_function, reconstruction_, param);
    range_vector.set_all(0.);
    advection_operator_.apply(reconstruction_.dofs().vector, range_vector);
  }

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override final
  {
    return reconstruction_operator_.source_space();
  }

  const AssemblyGridViewType& assembly_grid_view() const override final
  {
    return advection_operator_.assembly_grid_view();
  }

  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    reconstruction_operator_.apply(source_vector, reconstruction_, param);
    range_vector.set_all(0.);
    advection_operator_.apply(reconstruction_.dofs().vector, range_vector);
  }

  /// \}

protected:
  const AdvectionOperatorType& advection_operator_;
  const ReconstructionOperatorType& reconstruction_operator_;
  mutable VectorType reconstruction_;
}; // class AdvectionWithReconstructionOperator<...>


template <class AdvectionOperatorImp, class ReconstructionOperatorImp>
class AdvectionWithPointwiseReconstructionOperator
  : public OperatorInterface<typename AdvectionOperatorImp::AGV,
                             AdvectionOperatorImp::s_r,
                             AdvectionOperatorImp::s_rC,
                             AdvectionOperatorImp::r_r,
                             AdvectionOperatorImp::r_rC,
                             typename AdvectionOperatorImp::F,
                             typename AdvectionOperatorImp::M,
                             typename AdvectionOperatorImp::SGV,
                             typename AdvectionOperatorImp::RGV>
{
public:
  using ThisType = AdvectionWithPointwiseReconstructionOperator;
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::AGV,
                                     AdvectionOperatorImp::s_r,
                                     AdvectionOperatorImp::s_rC,
                                     AdvectionOperatorImp::r_r,
                                     AdvectionOperatorImp::r_rC,
                                     typename AdvectionOperatorImp::F,
                                     typename AdvectionOperatorImp::M,
                                     typename AdvectionOperatorImp::SGV,
                                     typename AdvectionOperatorImp::RGV>;

  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AdvectionOperatorType = AdvectionOperatorImp;
  using ReconstructionOperatorType = ReconstructionOperatorImp;

  AdvectionWithPointwiseReconstructionOperator(const AdvectionOperatorType& advection_operator,
                                               const ReconstructionOperatorType& reconstruction_operator)
    : BaseType(advection_operator.parameter_type() + reconstruction_operator.parameter_type())
    , advection_operator_(advection_operator)
    , reconstruction_operator_(reconstruction_operator)
    , reconstructed_values_(source_space().grid_view().indexSet().size(0))
    , reconstructed_function_(source_space().grid_view(), reconstructed_values_)
  {}

  // pull in methods from various base classes
  using BaseType::apply;

  /// \name Required by ForwardOperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override final
  {
    return advection_operator_.range_space();
  }

  bool linear() const override final
  {
    return false;
  }

  // avoid non-optimal default implementation in OperatorInterface
  void apply(SourceFunctionType source_function,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_range(range_vector);
    reconstruction_operator_.apply(source_function, reconstructed_function_, param);
    range_vector.set_all(0.);
    advection_operator_.apply(reconstructed_function_, range_vector, param);
  }

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override final
  {
    return advection_operator_.source_space();
  }

  const AssemblyGridViewType& assembly_grid_view() const override final
  {
    return advection_operator_.assembly_grid_view();
  }

  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    reconstruction_operator_.apply(source_vector, reconstructed_function_, param);
    range_vector.set_all(0.);
    advection_operator_.apply(reconstructed_function_, range_vector, param);
  }

  /// \}

  const AdvectionOperatorType& advection_operator_;
  const ReconstructionOperatorType& reconstruction_operator_;
  typename ReconstructionOperatorType::ReconstructedValuesType reconstructed_values_;
  mutable typename ReconstructionOperatorType::ReconstructedFunctionType reconstructed_function_;
}; // class AdvectionWithReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_WITH_RECONSTRUCTION_HH
