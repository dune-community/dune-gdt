#ifndef DUNE_DIRICHLETCONSTRAINTS_HH
#define DUNE_DIRICHLETCONSTRAINTS_HH

#include <dune/fem/function/common/scalarproducts.hh>
#include <string>
#include <sstream>

namespace Dune {

// TODO template parameters should be replaced by FunctionalType
template <class DiscreteFunctionSpace>
class DirichletConstraints
{

public:
  typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

  //! type of grid partition
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  //! type of grid
  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  // types for boundary treatment
  // ----------------------------
  typedef typename DiscreteFunctionSpaceType::MapperType MapperType;
  typedef SlaveDofs<DiscreteFunctionSpaceType, MapperType> SlaveDofsType;
  typedef typename SlaveDofsType::SingletonKey SlaveDofsKeyType; /*@\label{poi:singleton}@*/
  typedef SingletonList<SlaveDofsKeyType, SlaveDofsType> SlaveDofsProviderType;

  // TODO change constructor, remove problem
  DirichletConstraints(const DiscreteFunctionSpaceType& space)
    : space_(space)
    , slaveDofs_(getSlaveDofs(space_))
  {
  }

  template <class LinearOperator>
  void applyToOperator(LinearOperator& linearOperator) const
  {
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    typedef typename IteratorType::Entity EntityType;

    const IteratorType end = space_.end();
    for (IteratorType it = space_.begin(); it != end; ++it) {
      const EntityType& entity = *it;
      // if entity has boundary intersections
      if (entity.hasBoundaryIntersections()) {
        applyToOperatorLocal(linearOperator, entity);
      }
    }
  }


  /*! treatment of Dirichlet-DoFs for one entity
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[in]  entity  entity to perform Dirichlet treatment on
   */
  template <class LinearOperator, class EntityType> /*@LST0S@*/
  void applyToOperatorLocal(LinearOperator& linearOperator, const EntityType& entity) const
  { /*@LST0E@*/
    const int faceCodim = 1;
    typedef typename DiscreteFunctionSpaceType::LagrangePointSetType LagrangePointSetType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename LagrangePointSetType::template Codim<faceCodim>::SubEntityIteratorType FaceDofIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename LinearOperator::LocalMatrixType LocalMatrixType;

    const DiscreteFunctionSpaceType& dfSpace = space_;
    const GridPartType& gridPart             = dfSpace.gridPart();

    const LagrangePointSetType& lagrangePointSet /*@\label{poi:lagPoSet}@*/
        = dfSpace.lagrangePointSet(entity);

    // get slave dof structure (for parallel runs)   /*@LST0S@*/
    SlaveDofsType& slaveDofs = this->slaveDofs();
    const int numSlaveDofs   = slaveDofs.size();

    // get local matrix from matrix object
    LocalMatrixType localMatrix = linearOperator.localMatrix(entity, entity);

    // loop over all intersections and find dirichlet
    // boundaries
    const IntersectionIteratorType endit = gridPart.iend(entity);
    for (IntersectionIteratorType it = gridPart.ibegin(entity); it != endit; ++it) { /*@LST0E@*/
      // get intersection
      const IntersectionType& intersection = *it;

      // skip non-boundary elements
      if (!intersection.boundary()) /*@\label{poi:skipNonBd}@*/
        continue;

      // get local face number of boundary intersection
      const int face = intersection.indexInInside();

      // get iterator over all local dofs on this face
      FaceDofIteratorType faceIt          = lagrangePointSet.template beginSubEntity<faceCodim>(face);
      const FaceDofIteratorType faceEndIt = lagrangePointSet.template endSubEntity<faceCodim>(face);

      // iterate over face dofs and set unit row    /*@LST0S@*/
      for (; faceIt != faceEndIt; ++faceIt) {
        const int localDof = *faceIt;
        // clear the whole row
        localMatrix.clearRow(localDof);

        // set diagonal to 1
        localMatrix.set(localDof, localDof, 1);

        // clear all other columns
        const int globalDof = dfSpace.mapToGlobal(entity, localDof); /*@\label{poi:clslv0}@*/

        // cancel all slave dofs (for parallel runs)
        for (int i = 0; i < numSlaveDofs; ++i) {
          // slave dofs are canceled
          if (globalDof == slaveDofs[i])
            localMatrix.set(localDof, localDof, 0);
        } /*@\label{poi:clslv1}@*/
      }
    }
  }


public:
#if 0
    //! set the dirichlet points to exact values
    template< class EntityType, class DiscreteFunctionType >
    bool boundaryTreatment( const EntityType &entity,
                            const ProblemType& problem,
                            DiscreteFunctionType &rhs,
                            DiscreteFunctionType &solution ) const 
    {                                                                 /*@LST0E@*/
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
        DiscreteSpaceType;
      typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

      typedef typename DiscreteSpaceType :: LagrangePointSetType
        LagrangePointSetType;
      typedef typename DiscreteSpaceType :: GridPartType GridPartType;

      const int faceCodim = 1;
      typedef typename GridPartType :: IntersectionIteratorType
        IntersectionIteratorType;
      typedef typename LagrangePointSetType
        :: template Codim< faceCodim > :: SubEntityIteratorType
        FaceDofIteratorType;

      const DiscreteSpaceType &dfSpace = rhs.space();
      const GridPartType &gridPart = dfSpace.gridPart();

      // get local functions of data 
      LocalFunctionType rhsLocal = rhs.localFunction( entity );
      LocalFunctionType solutionLocal = solution.localFunction( entity );

      bool hasDirichletBoundary = false;

      typedef typename EntityType :: Geometry Geometry; 
      const Geometry& geo = entity.geometry();

      const LagrangePointSetType &lagrangePointSet = dfSpace.lagrangePointSet( entity );

      IntersectionIteratorType it = gridPart.ibegin( entity );
      const IntersectionIteratorType endit = gridPart.iend( entity );
      for( ; it != endit; ++it )
      {
        typedef typename IntersectionIteratorType :: Intersection IntersectionType;
        const IntersectionType& intersection = *it;

        // if intersection is with boundary, adjust data  
        if( intersection.boundary() )
        {
          hasDirichletBoundary = true;

          // get face number of boundary intersection 
          const int face = intersection.indexInInside();
          // get dof iterators 
          FaceDofIteratorType faceIt
            = lagrangePointSet.template beginSubEntity< faceCodim >( face );
          const FaceDofIteratorType faceEndIt
            = lagrangePointSet.template endSubEntity< faceCodim >( face );
          for( ; faceIt != faceEndIt; ++faceIt )
          {
            // get local dof number 
            const int localDof = *faceIt;

            typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
            // get global coordinate of point on boundary 
            const DomainType global = geo.global( lagrangePointSet.point( localDof ) );

            // check whether Dirichlet boundary or not 
            if( ! problem.dirichletBoundary(intersection.boundaryId(), global ) )
            {
              continue;
            }

            // evaluate boundary data
            typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
            RangeType phi;
            problem.g( global, phi );

            // adjust right hand side and solution data 
            rhsLocal[ localDof ] = phi[ 0 ];
            solutionLocal[ localDof ] = phi[ 0 ];
          }
        }
      }

      return hasDirichletBoundary;
    }                                                                 /*@LST0S@*/

#endif


protected:
  //! pointer to slave dofs
  const DiscreteFunctionSpaceType& space_;

  SlaveDofsType* const slaveDofs_;

  // return slave dofs                          /*@\label{poi:slavedof0}@*/
  static SlaveDofsType* getSlaveDofs(const DiscreteFunctionSpaceType& space)
  {
    SlaveDofsKeyType key(space, space.mapper());
    return &(SlaveDofsProviderType::getObject(key));
  }

  // return reference to slave dofs
  SlaveDofsType& slaveDofs() const
  {
    slaveDofs_->rebuild();
    return *slaveDofs_;
  } /*@\label{poi:slavedof1}@*/ /*@LST0E@*/
};

} // end namespace Dune
#endif
