#ifndef DIRICHLET_SHETXCOC
#define DIRICHLET_SHETXCOC

namespace Dune {
namespace Fem {
namespace Functional {
namespace Constraints {

template <class DiscFuncSpace>
class Dirichlet
{
public:
  typedef DiscFuncSpace DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType GridPartType;

  static const int griddim = GridPartType::GridType::dimension;

  typedef Constraints::LocalDefault<double, griddim, griddim> LocalConstraintsType;

public:
  template <class DomainConstraint>
  Dirichlet(DiscFuncSpace& space, DomainConstraint& domainConstraint)
    : space_(space)
    , gridPart_(space.gridPart())
  {
  }


  template <class Entity>
  const LocalConstraintsType local(const Entity& en)
  {
    typedef typename DiscFuncSpace::BaseFunctionSet BFS;

    const BFS& bfs    = space_.baseFunctionSet(en);
    const int numCols = bfs.numBaseFunctions();
    LocalConstraintsType lc(numCols);

    typedef typename DFS::Iterator ItType;
    typedef typename Entity::LeafIntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;

    int numRows = 0;
    ItType it = space_.begin();
    for (; it != space_.end(); it++) {
      const Entity& en = *it;

      IntersectionIterator iit = en.ileafbegin();
      for (; iit != en.ileafend(); iit++) {
        const Intersection& ii = *iit;
        if (ii.boundary()) {
          lc.setRowDofs(numRows, space_.mapToGlobal(en, numRows));
          for (unsigned int i = 0; i < numCols; i++) {
            lc.setColumnDofs(i, space_.mapToGlobal(en, i));
            lc.setLocalMatrix(numRows, i) = 0;
          }
          lc.setLocalMatrix(numRows, numRows) = 1;
          ++numRows;
        }
      }
    }
    lc.setRowDofsSize(numRows);
    return lc;
  }

private:
  DiscreteFunctionSpace& space_;
  GridPartType& gridPart_;
}; // end class Dirichlet

} // end of namespace Constraints
} // end of namespace Functional
} // end of namespace Fem
} // end of namespace Dune

#endif /* end of include guard: DIRICHLET_SHETXCOC */
