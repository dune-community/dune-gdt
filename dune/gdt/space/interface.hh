// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_INTERFACE_HH
#define DUNE_GDT_SPACE_INTERFACE_HH

#include <dune/common/dynvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/la/container/pattern.hh>

#include "constraints.hh"
#include "../discretefunction/visualize.hh"

// needs to be here, since one of the above includes the config.h
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune {
namespace GDT {


template< class Traits >
class SpaceInterface
{
public:
  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::GridPartType   GridPartType;
  static const int                        polOrder = Traits::polOrder;
  typedef typename GridPartType::ctype    DomainFieldType;
  static const unsigned int               dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::EntityType           EntityType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  const GridPartType& gridPart() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
    return asImp().gridPart();
  }

  const BackendType& backend() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().backend());
    return asImp().backend();
  }

  bool continuous() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().continuous());
    return asImp().continuous();
  }

  const MapperType& mapper() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().mapper());
    return asImp().mapper();
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().baseFunctionSet(entity));
    return asImp().baseFunctionSet(entity);
  }

  template< class ConstraintsType >
  void localConstraints(const EntityType& entity, ConstraintsType& ret) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().localConstraints(entity, ret));
    asImp().localConstraints(entity, ret);
  }

  PatternType* computePattern() const
  {
    return computePattern(*this);
  }

  template< class OtherSpaceType >
  PatternType* computePattern(const OtherSpaceType& other) const
  {
    // check type of the other space
    return computePattern(gridPart(), other);
  }

  /**
   *  \brief  computes a sparsity pattern, where this space is the test space (rows/outer) and the other space is the
   *          ansatz space (cols/inner)
   */
  template< class LocalGridPartType, class O >
  PatternType* computePattern(const LocalGridPartType& localGridPart,
                              const SpaceInterface< O >& otherSpace) const
  {
    PatternType* ret = new PatternType(mapper().size());
    PatternType& pattern = *ret;
    // walk the grid part
    for (typename LocalGridPartType::template Codim< 0 >::IteratorType entityIt = localGridPart.template begin< 0 >();
         entityIt != localGridPart.template end< 0 >();
         ++entityIt) {
      const typename LocalGridPartType::template Codim< 0 >::EntityType& entity = *entityIt;
      // get basefunctionsets
      const auto testBase = baseFunctionSet(entity);
      const auto ansatzBase = otherSpace.baseFunctionSet(entity);
      Dune::DynamicVector< size_t > globalRows(testBase.size(), 0);
      mapper().globalIndices(entity, globalRows);
      Dune::DynamicVector< size_t > globalCols(ansatzBase.size(), 0);
      otherSpace.mapper().globalIndices(entity, globalCols);
      for (size_t ii = 0; ii < testBase.size(); ++ii) {
        auto& columns = pattern.inner(globalRows[ii]);
        for (size_t jj = 0; jj < ansatzBase.size(); ++jj) {
          columns.insert(globalCols[jj]);
        }
      }
    } // walk the grid part
    return ret;
  } // ... computePattern(...)

  template< class VectorType >
  void visualize(const VectorType& vector, const std::string filename, const std::string name = "vector") const
  {
    typedef DiscreteFunction::VisualizationAdapter< derived_type, VectorType > DiscreteFunctionType;
    const auto function = std::make_shared< const DiscreteFunctionType >(asImp(), vector, name);
    typedef typename Dune::VTKWriter< typename GridPartType::GridViewType > VTKWriterType;
    VTKWriterType vtk_writer(gridPart().gridView(), Dune::VTK::nonconforming);
    vtk_writer.addVertexData(function);
    vtk_writer.write(filename);
  }

  void visualize(const std::string filename_prefix = "") const
  {
    std::string filename = filename_prefix;
    if (filename_prefix.empty()) {
      filename = "raviartthomasspace.femlocalfunctions";
//                 + Stuff::Common::toString(dimDomain) + "_to_" + Stuff::Common::toString(dimRange);
      if (dimRangeCols > 1)
        filename += "x" + Stuff::Common::toString(dimRangeCols);
    }
    typedef typename Dune::VTKWriter< typename GridPartType::GridViewType > VTKWriterType;
    VTKWriterType vtk_writer(gridPart().gridView(), Dune::VTK::nonconforming);
    for (size_t ii = 0; ii < mapper().maxNumDofs(); ++ii) {
      typedef DiscreteFunction::BasisVisualization< typename Traits::derived_type > BasisFunctionType;
      const auto iith_baseFunction = std::make_shared< BasisFunctionType >(asImp(),
                                                                           ii,
                                                                           Stuff::Common::toString(ii) + "th basis");
      vtk_writer.addVertexData(iith_baseFunction);
    }
    vtk_writer.write(filename);
  }

  derived_type& asImp()
  {
    return static_cast< derived_type& >(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class SpaceInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_INTERFACE_HH
