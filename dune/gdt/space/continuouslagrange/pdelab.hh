// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_CONTINUOUSLAGRANGE_PDELAB_HH
#define DUNE_GDT_SPACE_CONTINUOUSLAGRANGE_PDELAB_HH

#include <memory>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

#include <dune/grid/common/capabilities.hh>

#ifdef HAVE_DUNE_PDELAB
# include <dune/pdelab/finiteelementmap/pkfem.hh>
# include <dune/pdelab/finiteelementmap/qkfem.hh>
# include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#endif // HAVE_DUNE_PDELAB

#include "../../mapper/pdelab.hh"
#include "../../basefunctionset/pdelab.hh"

#include "../interface.hh"
#include "../constraints.hh"

namespace Dune {
namespace GDT {
namespace ContinuousLagrangeSpace {

#ifdef HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabWrapper
{
  static_assert((Dune::AlwaysFalse< GridPartImp >::value), "Untested for this combination of dimensions!");
};


/**
 *  \brief Traits class for ContinuousLagrangeSpace::PdelabWrapper.
 */
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabWrapperTraits
{
public:
  typedef PdelabWrapper< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  typedef GridPartImp                   GridPartType;
  static const int                      polOrder = polynomialOrder;
  static_assert(polOrder >= 1, "Wrong polOrder given!");
private:
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridPartType::dimension;
public:
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = rangeDim;
  static const unsigned int             dimRangeCols = rangeDimCols;
private:
  template< class G, bool single_geom, bool is_simplex, bool is_cube >
  struct FeMap
  {
    static_assert(Dune::AlwaysFalse< G >::value,
                  "This space is only implemented for either fully simplicial or fully cubic grids!");
  };
  template< class G >
  struct FeMap< G, true, true, false >
  {
    typedef PDELab::PkLocalFiniteElementMap
      < typename GridPartType::GridViewType, DomainFieldType, RangeFieldType, polOrder> Type;
  };
  template< class G >
  struct FeMap< G, true, false, true >
  {
    typedef PDELab::QkLocalFiniteElementMap
      < typename GridPartType::GridViewType, DomainFieldType, RangeFieldType, polOrder> Type;
  };
  typedef typename GridPartType::GridType GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType< GridType >::v;
  static const bool simplicial_ = (Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
                                   == GenericGeometry::SimplexTopology< dimDomain >::type::id);
  static const bool cubic_ = (Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
                              == GenericGeometry::CubeTopology< dimDomain >::type::id);
  typedef typename FeMap< GridType, single_geom_, simplicial_, cubic_ >::Type FEMapType;
public:
  typedef PDELab::GridFunctionSpace< typename GridPartType::GridViewType, FEMapType > BackendType;
  typedef Mapper::PdelabPkQk< BackendType > MapperType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef BaseFunctionSet::PdelabWrapper
      < BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols >
    BaseFunctionSetType;
private:
  friend class PdelabWrapper< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >;
}; // class SpaceWrappedFemContinuousLagrangeTraits


template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class PdelabWrapper< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
  : public SpaceInterface< PdelabWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > >
{
  typedef SpaceInterface< PdelabWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > > BaseType;
  typedef PdelabWrapper< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >                         ThisType;
public:
  typedef PdelabWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > Traits;

  typedef typename Traits::GridPartType GridPartType;
  static const int                      polOrder = Traits::polOrder;
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridPartType::dimension;
  typedef FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::EntityType           EntityType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  PdelabWrapper(const std::shared_ptr< const GridPartType >& gridP)
    : gridPart_(gridP)
    , fe_map_(std::make_shared< FEMapType >(gridPart_->gridView()))
    , backend_(std::make_shared< BackendType >(const_cast< GridPartType& >(*gridPart_).gridView(),
                                               *fe_map_))
    , mapper_(std::make_shared< MapperType >(*backend_))
    , tmpMappedRows_(mapper_->maxNumDofs())
    , tmpMappedCols_(mapper_->maxNumDofs())
  {}

  PdelabWrapper(const ThisType& other)
    : gridPart_(other.gridPart_)
    , fe_map_(other.fe_map_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
    , tmpMappedRows_(mapper_->maxNumDofs())
    , tmpMappedCols_(mapper_->maxNumDofs())
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridPart_ = other.gridPart_;
      fe_map_ = other.fe_map_;
      backend_ = other.backend_;
      mapper_ = other.mapper_;
      tmpMappedRows_.resize(mapper_->maxNumDofs());
      tmpMappedCols_.resize(mapper_->maxNumDofs());
    }
    return *this;
  }

  ~PdelabWrapper() {}

  const std::shared_ptr< const GridPartType >& gridPart() const
  {
    return gridPart_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  bool continuous() const
  {
    return true;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return BaseFunctionSetType(*backend_, entity, polOrder);
  }

  template< class R >
  void localConstraints(const EntityType& /*entity*/, Constraints::LocalDefault< R >& /*ret*/) const
  {
    static_assert((Dune::AlwaysFalse< R >::value), "Not implemented for arbitrary constraints!");
  }

  void localConstraints(const EntityType& entity,
                        Constraints::Dirichlet
                          < typename GridPartType::IntersectionType, RangeFieldType, true >& ret) const
  {
    static_assert(dimDomain == 2, "Not tested for other dimensions!");
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    const std::set< size_t > localDirichletDofs = this->findLocalDirichletDoFs(entity, ret.gridBoundary());
    const size_t numRows = localDirichletDofs.size();
    if (numRows > 0) {
      const size_t numCols = mapper_->numDofs(entity);
      ret.setSize(numRows, numCols);
      mapper_->globalIndices(entity, tmpMappedRows_);
      mapper_->globalIndices(entity, tmpMappedCols_);
      size_t localRow = 0;
      const RangeFieldType zero(0);
      const RangeFieldType one(1);
      for (auto localDirichletDofIt = localDirichletDofs.begin();
           localDirichletDofIt != localDirichletDofs.end();
           ++localDirichletDofIt) {
        const size_t& localDirichletDofIndex = * localDirichletDofIt;
        ret.globalRow(localRow) = tmpMappedRows_[localDirichletDofIndex];
        for (size_t jj = 0; jj < ret.cols(); ++jj) {
          ret.globalCol(jj) = tmpMappedCols_[jj];
          if (tmpMappedCols_[jj] == tmpMappedRows_[localDirichletDofIndex])
            ret.value(localRow, jj) = one;
          else
            ret.value(localRow, jj) = zero;
        }
        ++localRow;
      }
    } else {
      ret.setSize(0, 0);
    }
  } // ... localConstraints(..., Dirichlet< ..., true >)

  void localConstraints(const EntityType& entity,
                        Constraints::Dirichlet
                          < typename GridPartType::IntersectionType, RangeFieldType, false >& ret) const
  {
    static_assert(dimDomain == 2, "Not tested for other dimensions!");
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    const std::set< size_t > localDirichletDofs = this->findLocalDirichletDoFs(entity, ret.gridBoundary());
    const size_t numRows = localDirichletDofs.size();
    if (numRows > 0) {
      const size_t numCols = mapper_->numDofs(entity);
      ret.setSize(numRows, numCols);
      mapper_->globalIndices(entity, tmpMappedRows_);
      mapper_->globalIndices(entity, tmpMappedCols_);
      size_t localRow = 0;
      const RangeFieldType zero(0);
      for (auto localDirichletDofIt = localDirichletDofs.begin();
           localDirichletDofIt != localDirichletDofs.end();
           ++localDirichletDofIt) {
        const size_t& localDirichletDofIndex = * localDirichletDofIt;
        ret.globalRow(localRow) = tmpMappedRows_[localDirichletDofIndex];
        for (size_t jj = 0; jj < ret.cols(); ++jj) {
          ret.globalCol(jj) = tmpMappedCols_[jj];
          ret.value(localRow, jj) = zero;
        }
        ++localRow;
      }
    } else {
      ret.setSize(0, 0);
    }
  } // ... localConstraints(..., Dirichlet< ..., false >)

  using BaseType::computePattern;

  template< class LocalGridPartType, class OtherSpaceType >
  PatternType* computePattern(const LocalGridPartType& localGridPart, const OtherSpaceType& otherSpace) const
  {
    return BaseType::computeCodim0Pattern(localGridPart, otherSpace);
  }

  std::vector< DomainType > lagrange_points(const EntityType& entity) const
  {
    // check
    static_assert(polOrder == 1, "Not yet implemented for other polynomial orders!");
    // get the basis and reference element
    const auto basis = baseFunctionSet(entity);
    const auto& reference_element = ReferenceElements< DomainFieldType, dimDomain >::general(entity.type());
    const int num_vertices = reference_element.size(dimDomain);
    assert(num_vertices >= 0);
    assert(size_t(num_vertices) == basis.size() && "This should not happen with polOrder 1!");
    // prepare return vector
    std::vector< DomainType > local_vertices(num_vertices, DomainType(0));
    if (this->tmp_basis_values_.size() < basis.size())
      this->tmp_basis_values_.resize(basis.size());
    // loop over all vertices
    for (int ii = 0; ii < num_vertices; ++ii) {
      // get the local coordinate of the iith vertex
      const auto local_vertex = reference_element.position(ii, dimDomain);
      // evaluate the basefunctionset
      basis.evaluate(local_vertex, this->tmp_basis_values_);
      // find the basis function that evaluates to one here (has to be only one!)
      size_t found = 0;
      for (size_t jj = 0; jj < basis.size(); ++jj)
        if (Dune::Stuff::Common::FloatCmp::eq(this->tmp_basis_values_[jj][0], RangeFieldType(1))) {
          ++found;
          local_vertices[jj] = local_vertex;
        }
      assert(found == 1 && "This must not happen for polOrder 1!");
    }
    return local_vertices;
  } // ... lagrange_points(...)

private:
  std::shared_ptr< const GridPartType > gridPart_;
  std::shared_ptr< const FEMapType > fe_map_;
  std::shared_ptr< const BackendType > backend_;
  std::shared_ptr< const MapperType > mapper_;
  mutable Dune::DynamicVector< size_t > tmpMappedRows_;
  mutable Dune::DynamicVector< size_t > tmpMappedCols_;
}; // class PdelabWrapper< ..., 1 >


#else // HAVE_DUNE_PDELAB


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabWrapper
{
  static_assert((Dune::AlwaysFalse< GridPartImp >::value), "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace ContinuousLagrangeSpace
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_CONTINUOUSLAGRANGE_PDELAB_HH
