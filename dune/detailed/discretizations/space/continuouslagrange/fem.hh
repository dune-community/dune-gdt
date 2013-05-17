#ifndef DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONTINUOUS_LAGRANGE_FEM_HH
#define DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONTINUOUS_LAGRANGE_FEM_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrangespace.hh>

#include "../../mapper/fem.hh"
#include "../../basefunctionset/fem.hh"

#include "../interface.hh"
#include "../constraints.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace ContinuousLagrangeSpace {


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDimRows = 1, int rangeDimCols = 1>
class FemWrapper;


// forward, to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDimRows = 1, int rangeDimCols = 1>
class FemWrapperTraits;


/**
 *  \brief Traits class for ContinuousLagrangeSpace for dimRange 1x1.
 */
template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class FemWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
{
public:
  typedef FemWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> derived_type;
  typedef GridPartImp GridPartType;
  static const int polOrder = polynomialOrder;
  dune_static_assert((polOrder >= 1), "ERROR: wrong polOrder given!");

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = 1;
  static const unsigned int dimRangeCols = 1;

private:
  typedef Dune::Fem::FunctionSpace<DomainFieldType, RangeFieldType, dimDomain, dimRangeRows> FunctionSpaceType;

public:
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType, polOrder> BackendType;
  typedef Mapper::FemDofWrapper<typename BackendType::MapperType> MapperType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef BaseFunctionSet::FemWrapper<typename BackendType::BaseFunctionSetType, EntityType, DomainFieldType, dimDomain,
                                      RangeFieldType, dimRange, dimRangeCols> BaseFunctionSetType;
}; // class SpaceWrappedFemContinuousLagrangeTraits< ..., 1, 1 >


template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class FemWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public SpaceInterface<FemWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>>
{
public:
  typedef FemWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::GridPartType GridPartType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const int polOrder           = Traits::polOrder;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRangeRows = Traits::dimRangeRows;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

public:
  FemWrapper(const GridPartType& gridP)
    : gridPart_(gridP)
    , backend_(const_cast<GridPartType&>(gridPart_))
    , mapper_(backend_.mapper())
    , tmpMappedRows_(mapper_.maxNumDofs())
    , tmpMappedCols_(mapper_.maxNumDofs())
  {
  }

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  bool continuous() const
  {
    return true;
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return BaseFunctionSetType(backend_, entity);
  }

  template <class R>
  void localConstraints(const EntityType& /*entity*/, Constraints::LocalDefault<R>& /*ret*/) const
  {
    dune_static_assert(Dune::AlwaysFalse<R>::value, "ERROR: not implemented for arbitrary constraints!");
  }

  void localConstraints(const EntityType& entity,
                        Constraints::Dirichlet<typename GridPartType::GridViewType, RangeFieldType>& ret) const
  {
    const auto& lagrangePointSet = backend_.lagrangePointSet(entity);
    std::set<size_t> localDirichletDofs;
    // loop over all intersections
    const auto& gridBoundary     = ret.gridBoundary();
    const auto intersectionEndIt = gridPart_.iend(entity);
    for (auto intersectionIt = gridPart_.ibegin(entity); intersectionIt != intersectionEndIt; ++intersectionIt) {
      const auto& intersection = *intersectionIt;
      // only work on dirichlet intersections
      if (gridBoundary.dirichlet(intersection)) {
        // get local face number of boundary intersection
        const int intersectionIndex = intersection.indexInInside();
        // iterate over face dofs and set unit row
        const auto faceDofEndIt = lagrangePointSet.template endSubEntity<1>(intersectionIndex);
        for (auto faceDofIt = lagrangePointSet.template beginSubEntity<1>(intersectionIndex); faceDofIt != faceDofEndIt;
             ++faceDofIt) {
          const size_t localDofIndex = *faceDofIt;
          localDirichletDofs.insert(localDofIndex);
        } // iterate over face dofs and set unit row
      } // only work on dirichlet intersections
    } // loop over all intersections
    const BaseFunctionSetType basis = baseFunctionSet(entity);
    const size_t numRows = localDirichletDofs.size();
    if (numRows > 0) {
      const size_t numCols = basis.size();
      ret.setSize(numRows, numCols);
      mapper_.globalIndices(entity, tmpMappedRows_);
      mapper_.globalIndices(entity, tmpMappedCols_);
      size_t localRow = 0;
      const RangeFieldType zero(0);
      const RangeFieldType one(1);
      for (auto localDirichletDofIt = localDirichletDofs.begin(); localDirichletDofIt != localDirichletDofs.end();
           ++localDirichletDofIt) {
        const size_t& localDirichletDofIndex = *localDirichletDofIt;
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
  } // ... localConstraints(..., Dirichlet)

private:
  const GridPartType& gridPart_;
  const BackendType backend_;
  const MapperType mapper_;
  mutable Dune::DynamicVector<size_t> tmpMappedRows_;
  mutable Dune::DynamicVector<size_t> tmpMappedCols_;
}; // class FemWrapper< ..., 1, 1 >


} // namespace ContinuousLagrangeSpace
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_SPACE_CONTINUOUS_LAGRANGE_FEM_HH
