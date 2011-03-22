#ifndef DIRICHLET_SHETXCOC
#define DIRICHLET_SHETXCOC

namespace Dune {
namespace Functionals {
namespace Constraints {

/** @brief Constraints for Dirichlet values on the entire boundary domain
 *
 * This class implements constraints on the degrees of freedom on a @ref
 * Subspace::Linear "linear subspace" @endref.
 *
 * Constraints efficiently implement functionals @f$ f_i:\cal X_H \to
 * \mathbb{R}@f$ for each degree of freedom @f$i\in C \subset
 * \{1,\dots,H\}@f$, where $H$ is the number of degree of freedoms in the
 * underlying discrete function space @f$\cal X_H\@f$. In case of this Dirichlet constraints class, the set
 * @f$\{1,\dots,C\}@f$ includes all
 *
 * @note The Dirichlet constraints only make sense on a finite element space,
 * not on a discontinuous discrete function space.
 *
 * @tparam DiscFuncSpace discrete function space on which constraints shall
 * be applied.
 */
template <class DiscFuncSpace>
class Dirichlet
{
public:
  //! Discrete function space on which the Dirichlet constraints are applied
  typedef DiscFuncSpace DiscreteFunctionSpace;

  //! Underlying grid part
  typedef typename DiscreteFunctionSpace::GridPartType GridPartType;

  //! Dimension of the grid part
  static const int griddim = GridPartType::GridType::dimension;

  //! Return type of local() method, implementing the LocalConstraints
  //! interface
  typedef Constraints::LocalDefault<double, griddim, griddim> LocalConstraintsType;

public:
  /** Constructor for the Dirichlet constraints
   *
   *  @param space    discrete function space object on which the Dirichlet
   *                  constraints are applied
   */
  Dirichlet(DiscFuncSpace& space)
    : space_(space)
    , gridPart_(space.gridPart())
  {
  }


  /** return a local constraint object for the entity \a en
   *
   * @param en Entity for which the local constraints shall be compted
   *
   * @returns local constraints object (copyable).
   */
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
} // end of namespace Functionals
} // end of namespace Dune

#endif /* end of include guard: DIRICHLET_SHETXCOC */
