// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_INTERFACE_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_INTERFACE_HH

#include <memory>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/functions/default.hh>
#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace GDT {
namespace LinearElliptic {


template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class ProblemInterface
{
  typedef ProblemInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;

public:
  typedef EntityImp EntityType;
  typedef DomainFieldImp DomainFieldType;
  static const size_t dimDomain = domainDim;
  typedef RangeFieldImp RangeFieldType;
  static const size_t dimRange = rangeDim;

  typedef XT::Functions::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>
      DiffusionFactorType;
  typedef XT::Functions::
      LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain, dimDomain>
          DiffusionTensorType;
  typedef XT::Functions::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;

  virtual ~ProblemInterface()
  {
  }

  virtual const DiffusionFactorType& diffusion_factor() const = 0;

  virtual const DiffusionTensorType& diffusion_tensor() const = 0;

  virtual const FunctionType& force() const = 0;

  virtual const FunctionType& dirichlet() const = 0;

  virtual const FunctionType& neumann() const = 0;

  virtual const XT::Common::Configuration& grid_cfg() const = 0;

  virtual const XT::Common::Configuration& boundary_info_cfg() const = 0;

  template <class G>
  void visualize(const GridView<G>& grid_layer,
                 std::string filename,
                 const bool subsampling = true,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    auto vtk_writer =
        subsampling ? Dune::XT::Common::make_unique<SubsamplingVTKWriter<GridView<G>>>(grid_layer, VTK::nonconforming)
                    : Dune::XT::Common::make_unique<VTKWriter<GridView<G>>>(grid_layer, VTK::nonconforming);
    auto diffusion = XT::Functions::make_product(diffusion_factor(), diffusion_tensor(), "diffusion");
    add_function_visualization(grid_layer, diffusion_factor(), *vtk_writer);
    add_function_visualization(grid_layer, diffusion_tensor(), *vtk_writer);
    add_function_visualization(grid_layer, *diffusion, *vtk_writer);
    add_function_visualization(grid_layer, force(), *vtk_writer);
    add_function_visualization(grid_layer, dirichlet(), *vtk_writer);
    add_function_visualization(grid_layer, neumann(), *vtk_writer);
    vtk_writer->write(filename, vtk_output_type);
  } // ... visualize(...)

private:
  template <class GridLayerType, class F, class VTKWriterType>
  void
  add_function_visualization(const GridLayerType& /*grid_layer*/, const F& function, VTKWriterType& vtk_writer) const
  {
    typedef XT::Functions::VisualizationAdapterFunction<GridLayerType, F::dimRange, F::dimRangeCols>
        VisualizationAdapter;
    vtk_writer.addVertexData(std::make_shared<VisualizationAdapter>(function));
  }
}; // ProblemInterface


namespace internal {


template <class F>
struct is_problem_helper
{
  DXTC_has_typedef_initialize_once(EntityType);
  DXTC_has_typedef_initialize_once(DomainFieldType);
  DXTC_has_typedef_initialize_once(RangeFieldType);
  DXTC_has_static_member_initialize_once(dimDomain);
  DXTC_has_static_member_initialize_once(dimRange);

  static const bool is_candidate = DXTC_has_typedef(EntityType)<F>::value && DXTC_has_typedef(DomainFieldType)<F>::value
                                   && DXTC_has_typedef(RangeFieldType)<F>::value
                                   && DXTC_has_static_member(dimDomain)<F>::value
                                   && DXTC_has_static_member(dimRange)<F>::value;
}; // class is_problem_helper


} // namespace internal


template <class P, bool candidate = internal::is_problem_helper<P>::is_candidate>
struct is_problem : public std::is_base_of<ProblemInterface<typename P::EntityType,
                                                            typename P::DomainFieldType,
                                                            P::dimDomain,
                                                            typename P::RangeFieldType,
                                                            P::dimRange>,
                                           P>
{
};


template <class P>
struct is_problem<P, false> : public std::false_type
{
};


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_INTERFACE_HH
