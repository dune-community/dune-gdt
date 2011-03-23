#ifndef DIRICHLET_SHETXCOC
#define DIRICHLET_SHETXCOC

#include "localdefault.hh"

namespace Dune {
namespace Functionals {
namespace Constraints {

/** @brief Constraints for Dirichlet values on the entire boundary domain
 *
 * This class implements constraints on the degrees of freedom on a @link
 * Subspace::Linear "linear subspace" @endlink.
 *
 * Constraints efficiently implement functionals @f$ f_i:{\cal X}_H \to
 * \mathbb{R}@f$ for each degree of freedom @f$i\in {\cal C} \subset
 * \{1,\dots,H\}@f$, where @f$ H @f$ is the number of degree of freedoms in the
 * underlying discrete function space @f$ {\cal X}_H @f$. In case of this
 * Dirichlet constraints class, the set @f$ {\cal C} @f$ includes all degrees of
 * freedom "lying" on boundary edges and @f$ f_i(u) = 0 @f$ for all @f$ i=\cal
 * C @f$.
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

  //! @brief dimension of the grid part
  static const int griddim = GridPartType::GridType::dimension;

  //! @brief Return type of local() method, implementing the LocalConstraints
  //! interface
  typedef Constraints::LocalDefault<double, griddim, griddim> LocalConstraintsType;

public:
  /** @brief Constructor for the Dirichlet constraints
   *
   *  @param space    discrete function space object on which the Dirichlet
   *                  constraints are applied
   */
  Dirichlet(DiscFuncSpace& space)
    : space_(space)
    , gridPart_(space.gridPart())
  {
  }


  /** @brief returns a local constraint object for the entity \a en
   *
   * @param en Entity for which the local constraints shall be compted
   *
   * @returns local constraints object (copyable). c.f.
   * Constraints::LocalDefault for more details on the return type
   * implementation.
   */
  template <class Entity>
  const LocalConstraintsType local(const Entity& en)
  {
    typedef typename DiscFuncSpace::BaseFunctionSetType BFS;

    const BFS& bfs    = space_.baseFunctionSet(en);
    const int numCols = bfs.numBaseFunctions();
    LocalConstraintsType lc(numCols);

    typedef typename DiscFuncSpace::IteratorType ItType;
    typedef typename Entity::LeafIntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;

    unsigned int numRows = 0;
    ItType it = space_.begin();
    for (; it != space_.end(); ++it) {
      const Entity& en = *it;

      IntersectionIterator iit = en.ileafbegin();
      for (; iit != en.ileafend(); ++iit) {
        const Intersection& ii = *iit;
        if (ii.boundary()) {
          lc.setRowDofs(numRows, space_.mapToGlobal(en, numRows));
          for (unsigned int i = 0; i < numCols; ++i) {
            lc.setColumnDofs(i, space_.mapToGlobal(en, i));
            lc.setLocalMatrix(numRows, i, 0.0);
          }
          lc.setLocalMatrix(numRows, numRows, 1.0);
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
