// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_INTERFACE_HH
#define DUNE_GDT_SPACE_INTERFACE_HH

#include <memory>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/float_cmp.hh>

#include "constraints.hh"

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
  typedef typename GridPartType::ctype              DomainFieldType;
  static const unsigned int                         dimDomain = GridPartType::dimension;
  typedef FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::EntityType           EntityType;

  typedef typename GridPartType::IntersectionType           IntersectionType;
  typedef Stuff::GridboundaryInterface< IntersectionType >  BoundaryInfoType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  std::shared_ptr< const GridPartType > gridPart() const
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

  std::set< size_t > findLocalDirichletDoFs(const EntityType& entity, const BoundaryInfoType& boundaryInfo) const
  {
    static_assert(dimRange == 1, "Not implemented for this dimension!");
    static_assert(dimRangeCols == 1, "Not implemented for this dimension!");
    if (polOrder != 1) DUNE_THROW(Dune::NotImplemented, "This does not seem to work for higher orders!");
    // check
    assert(gridPart()->indexSet().contains(entity));
    // prepare
    std::set< size_t > localDirichletDofs;
    std::vector< typename BaseFunctionSetType::DomainType > dirichlet_vertices;
    // get all dirichlet vertices of this entity, therefore
    // * loop over all intersections
    const auto intersection_it_end = gridPart()->iend(entity);
    for (auto intersection_it = gridPart()->ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
      // only work on dirichlet ones
      const auto& intersection = *intersection_it;
      if (boundaryInfo.dirichlet(intersection)) {
        // and get the vertices of the intersection
        const auto geometry = intersection.geometry();
        for (size_t cc = 0; cc < geometry.corners(); ++cc)
          dirichlet_vertices.emplace_back(entity.geometry().local(geometry.corner(cc)));
      } // only work on dirichlet ones
    } // loop over all intersections

    // find the corresponding basis functions
    const auto basis = baseFunctionSet(entity);
    if (tmp_basis_values_.size() < basis.size())
      tmp_basis_values_.resize(basis.size());
    for (size_t cc = 0; cc < dirichlet_vertices.size(); ++cc) {
      basis.evaluate(dirichlet_vertices[cc], tmp_basis_values_);
      for (size_t jj = 0; jj < basis.size(); ++jj)
        if (Dune::Stuff::Common::FloatCmp::eq(tmp_basis_values_[jj][0], RangeFieldType(1)))
          localDirichletDofs.insert(jj);
    }

    // finish
    return localDirichletDofs;
  } // ... findLocalDirichletDoFs(...)

  template< class ConstraintsType >
  void localConstraints(const EntityType& entity, ConstraintsType& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().localConstraints(entity, ret));
    asImp().localConstraints(entity, ret);
  }

  PatternType* computePattern() const
  {
    return computePattern(*this);
  }

  template< class OtherSpaceType >
  PatternType* computePattern(const OtherSpaceType& other) const
  {
    return computePattern(*(gridPart()), other);
  }

  template< class LocalGridPartType, class OtherSpaceType >
  PatternType* computePattern(const LocalGridPartType& localGridPart,
                              const OtherSpaceType& otherSpace) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().computePattern(localGridPart, otherSpace));
    return asImp().computePattern(localGridPart, otherSpace);
  }

  /**
   *  \brief  computes a sparsity pattern, where this space is the test space (rows/outer) and the other space is the
   *          ansatz space (cols/inner)
   */
  template< class LocalGridPartType, class O >
  PatternType* computeCodim0Pattern(const LocalGridPartType& localGridPart,
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
  } // ... computeVolumePattern(...)

  /**
   *  \brief  computes a DG sparsity pattern, where this space is the test space (rows/outer) and the other space is the
   *          ansatz space (cols/inner)
   */
  template< class LocalGridPartType, class O >
  PatternType* computeCodim0AndCodim1Pattern(const LocalGridPartType& local_grid_part,
                                             const SpaceInterface< O >& other_space) const
  {
    // prepare
    PatternType* ret = new PatternType(mapper().size());
    PatternType& pattern = *ret;
    Dune::DynamicVector< size_t > global_rows(mapper().maxNumDofs(), 0);
    Dune::DynamicVector< size_t > global_cols(other_space.mapper().maxNumDofs(), 0);
    // walk the grid part
    const auto entity_it_end = local_grid_part.template end< 0 >();
    for (auto entity_it = local_grid_part.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get basefunctionsets
      const auto test_base_entity = baseFunctionSet(entity);
      const auto ansatz_base_entity = other_space.baseFunctionSet(entity);
      mapper().globalIndices(entity, global_rows);
      other_space.mapper().globalIndices(entity, global_cols);
      // compute entity/entity
      for (size_t ii = 0; ii < test_base_entity.size(); ++ii) {
        auto& columns = pattern.inner(global_rows[ii]);
        for (size_t jj = 0; jj < ansatz_base_entity.size(); ++jj) {
          columns.insert(global_cols[jj]);
        }
      }
      // walk the intersections
      const auto intersection_it_end = local_grid_part.iend(entity);
      for (auto intersection_it = local_grid_part.ibegin(entity);
           intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        // get the neighbour
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour = *neighbour_ptr;
          // get the basis
          const auto ansatz_base_neighbour = other_space.baseFunctionSet(neighbour);
          other_space.mapper().globalIndices(neighbour, global_cols);
          // compute entity/neighbour
          for (size_t ii = 0; ii < test_base_entity.size(); ++ii) {
            auto& columns = pattern.inner(global_rows[ii]);
            for (size_t jj = 0; jj < ansatz_base_neighbour.size(); ++jj) {
              columns.insert(global_cols[jj]);
            }
          }
        } // get the neighbour
      } // walk the intersections
    } // walk the grid part
    return ret;
  } // ... computeVolumeAndCouplingPattern(...)

  template< class LocalGridPartType, class O >
  PatternType* computeCodim1Pattern(const LocalGridPartType& local_grid_part,
                                    const SpaceInterface< O >& other_space) const
  {
    // prepare
    PatternType* ret = new PatternType(mapper().size());
    PatternType& pattern = *ret;
    Dune::DynamicVector< size_t > global_rows(mapper().maxNumDofs(), 0);
    Dune::DynamicVector< size_t > global_cols(other_space.mapper().maxNumDofs(), 0);
    // walk the grid part
    const auto entity_it_end = local_grid_part.template end< 0 >();
    for (auto entity_it = local_grid_part.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get basefunctionsets
      const auto test_base_entity = baseFunctionSet(entity);
      mapper().globalIndices(entity, global_rows);
      // walk the intersections
      const auto intersection_it_end = local_grid_part.iend(entity);
      for (auto intersection_it = local_grid_part.ibegin(entity);
           intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        // get the neighbour
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour = *neighbour_ptr;
          // get the basis
          const auto ansatz_base_neighbour = other_space.baseFunctionSet(neighbour);
          other_space.mapper().globalIndices(neighbour, global_cols);
          // compute entity/neighbour
          for (size_t ii = 0; ii < test_base_entity.size(); ++ii) {
            auto& columns = pattern.inner(global_rows[ii]);
            for (size_t jj = 0; jj < ansatz_base_neighbour.size(); ++jj) {
              columns.insert(global_cols[jj]);
            }
          }
        } // get the neighbour
      } // walk the intersections
    } // walk the grid part
    return ret;
  } // ... computeCodim1Pattern(...)

private:
  class BasisVisualization
    : public Dune::VTKFunction< typename GridPartType::GridViewType >
  {
    static_assert(dimRangeCols == 1, "Not implemented yet!");
  public:
    typedef typename BaseFunctionSetType::DomainType  DomainType;
    typedef typename BaseFunctionSetType::RangeType   RangeType;

    BasisVisualization(const derived_type& sp, const size_t ind, const std::string nm = "basis")
      : space_(sp)
      , index_(ind)
      , values_(space_.mapper().maxNumDofs(), RangeType(0))
      , name_(nm)
    {
      if (index_ >= space_.mapper().maxNumDofs())
        DUNE_THROW(Dune::RangeError,
                   "index has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << index_ << ")!");
    }

    virtual std::string name() const
    {
      return name_;
    }

    /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
    /* @{ */
    virtual int ncomps() const
    {
      return dimRange;
    }

    virtual double evaluate(int component, const EntityType& entity, const DomainType& xx) const
    {
      const auto baseFunctionSet = space_.baseFunctionSet(entity);
      if (component < 0)
        DUNE_THROW(Dune::RangeError, "component must not be negative (is " << component << ")!");
      if (component < int(baseFunctionSet.size())) {
        baseFunctionSet.evaluate(xx, values_);
        assert(component < int(values_.size()) && "This should not happen!");
        return values_[index_][component];
      } else if (component < int(space_.mapper().maxNumDofs()))
        return 0.0;
      else
        DUNE_THROW(Dune::RangeError,
                   "component has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << component << ")!");
    }
    /* @} */

  private:
    const derived_type& space_;
    const size_t index_;
    mutable std::vector< RangeType > values_;
    const std::string name_;
  }; // class BasisVisualization

public:
  void visualize(const std::string filename_prefix = "") const
  {
    std::string filename = filename_prefix;
    if (filename.empty()) {
      filename = "dune.gdt.space";
    }
    typedef typename Dune::VTKWriter< typename GridPartType::GridViewType > VTKWriterType;
    VTKWriterType vtk_writer(gridPart()->gridView(), Dune::VTK::nonconforming);
    for (size_t ii = 0; ii < mapper().maxNumDofs(); ++ii) {
      std::string number = "";
      if (ii == 1)
        number = "1st";
      else if (ii == 2)
        number = "2nd";
      else if (ii == 3)
        number = "3rd";
      else
        number = Stuff::Common::toString(ii) + "th";
      const auto iith_baseFunction = std::make_shared< BasisVisualization >(asImp(), ii, number + " basis");
      vtk_writer.addVertexData(iith_baseFunction);
    }
    vtk_writer.write(filename);
  } // ... visualize(...)

  derived_type& asImp()
  {
    return static_cast< derived_type& >(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }

private:
  mutable std::vector< typename BaseFunctionSetType::RangeType > tmp_basis_values_;
}; // class SpaceInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_INTERFACE_HH
