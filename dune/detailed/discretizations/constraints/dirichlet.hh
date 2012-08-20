#ifndef DUNE_DETAILED_DISCRETIZATIONS_CONSTRAINTS_DIRICHLET_HH
#define DUNE_DETAILED_DISCRETIZATIONS_CONSTRAINTS_DIRICHLET_HH

// system
#include <set>

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>

// local includes
#include "localdefault.hh"

// dune-stuff
#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace Constraints
{

/** @brief Constraints for zero Dirichlet values on the entire boundary domain
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
 * @ingroup Constraints
 *
 * @tparam DiscFuncSpace discrete function space on which constraints shall
 * be applied.
 */
template< class DiscreteFunctionSpaceImp, class BoundaryInfoImp = Dune::Stuff::Grid::BoundaryInfo::AllDirichlet >
class DirichletZero
{
public:
  //! Discrete function space on which the Dirichlet constraints are applied
  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef BoundaryInfoImp BoundaryInfoType;

  typedef DirichletZero< DiscreteFunctionSpaceType, BoundaryInfoType >
    ThisType;

  //! Underlying grid part
  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;

  //! @brief dimension of the grid part
  static const int griddim = GridPartType::GridType::dimension;

  //! @brief Return type of local() method, implementing the LocalConstraints
  //! interface
  typedef Dune::Detailed::Discretizations::Constraints::LocalDefault< double, 2*griddim, 6*griddim >
    LocalConstraintsType;

  /** @brief Constructor for the Dirichlet constraints
   *
   *  @param space    discrete function space object on which the Dirichlet
   *                  constraints are applied
   */
  DirichletZero(const DiscreteFunctionSpaceType& space,
                const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo = Dune::shared_ptr< const BoundaryInfoType >(new Dune::Stuff::Grid::BoundaryInfo::AllDirichlet()))
    : space_( space )
    , boundaryInfo_(boundaryInfo)
    , gridPart_( space.gridPart() )
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
  template< class Entity >
  LocalConstraintsType local( const Entity& entity ) const
  {
//    typedef typename DiscreteFunctionSpaceType::HostSpaceType
//      DiscFuncSpace;
//    typedef typename DiscreteFunctionSpaceType::LocalBaseFunctionSetType
//      LocalBaseFunctionSetType;
//    typedef typename DiscreteFunctionSpaceType::MapperType::LagrangePointSetType
//      LPSType;
    typedef LagrangePointSet< GridPartType, DiscreteFunctionSpaceType::polynomialOrder >
      LagrangePointSetType;
    typedef typename GridPartType::IntersectionIteratorType
      IntersectionIterator;
    typedef typename IntersectionIterator::Intersection
      Intersection;

    const int faceCodim = 1;
    typedef typename LagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType
       FaceDofIteratorType;

//    const LocalBaseFunctionSetType& localBaseFunctionSet = space_.localBaseFunctionSet( entity );
    const unsigned int numCols = space_.baseFunctionSet().local( entity ).size();
    LocalConstraintsType ret( numCols );

    const LagrangePointSetType& lagrangePointSet = space_.map().lagrangePointSet( entity );

/*    // get slave dof structure (for parallel runs)
 *    SlaveDofsType &slaveDofs = this->slaveDofs();
 *    const int numSlaveDofs = slaveDofs.size();*/

    // set of local boundary dofs
    std::set< unsigned int > localBoundaryDofs;

    // loop over all intersections
    const IntersectionIterator endIt = gridPart_.iend( entity );
    for( IntersectionIterator it = gridPart_.ibegin( entity ); it != endIt ; ++it )
    {
      // get intersection
      const Intersection& ii = *it;

      // only work on dirichlet intersections
      if (boundaryInfo_->dirichlet(ii)) {

        // get local face number of boundary intersection
        const int face = ii.indexInInside();

        // get iterator over all local dofs on this face
        FaceDofIteratorType faceIt
          = lagrangePointSet.template beginSubEntity< faceCodim >( face );
        const FaceDofIteratorType faceEndIt
          = lagrangePointSet.template endSubEntity< faceCodim >( face );

        // iterate over face dofs and set unit row
        for( ; faceIt != faceEndIt; ++faceIt )
        {
          const int localDof = *faceIt;

          localBoundaryDofs.insert( localDof );

  /*        // clear all other columns
   *        const int globalDof = dfSpace.mapToGlobal( entity, localDof );

   *        // cancel all slave dofs (for parallel runs)
   *        for( int i = 0; i < numSlaveDofs; ++i )
   *        {
   *          // slave dofs are canceled
   *          if( globalDof == slaveDofs[ i ] )
   *            localMatrix.set( localDof, localDof, 0 );
   *        }*/
        }
      } // only work on dirichlet intersections
    } // loop over all intersections


    /************************************************************************
     * iterate over local boundary dof set and fill local constraint matrix *
     ************************************************************************/
    typedef std::set< unsigned int >::const_iterator
      LBIterator;

    LBIterator lbend = localBoundaryDofs.end();
    unsigned int numRows = 0;
    for( LBIterator lbit = localBoundaryDofs.begin(); lbit != lbend; ++lbit, ++numRows )
    {
      const unsigned int localDof = *lbit;
      for( unsigned int i = 0; i < numCols; ++i )
      {
        if( numRows == 0 )
          ret.setColumnDofs( i, space_.map().toGlobal( entity, i ) );
        ret.setLocalMatrix( numRows, i, 0.0 );
      }
      ret.setLocalMatrix( numRows, localDof, 1.0 );
      ret.setRowDofs( numRows, space_.map().toGlobal( entity, localDof ) );
    }
    ret.setRowDofsSize( numRows );

    return ret;
  }

private:
  //! copy constructor
  DirichletZero( const ThisType& );

  //! assignment operator
  ThisType& operator=( const ThisType& );

  const DiscreteFunctionSpaceType& space_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const GridPartType& gridPart_;
}; // end class Dirichlet

} // end of namespace Constraints

} // namespace Discretizations

} // end of namespace Detailed

} // end of namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_CONSTRAINTS_DIRICHLET_HH
