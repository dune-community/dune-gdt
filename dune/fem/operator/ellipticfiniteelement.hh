#ifndef DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH
#define DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH

// dune fem includes
#include <dune/fem/function/common/function.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-fem-functionals includes
#include <dune/fem/common/localmatrix.hh>
#include <dune/fem/common/localvector.hh>
#include <dune/fem/functional/ltwo.hh>

namespace Dune {

namespace Functionals {

namespace Operator {

template <class DiscreteFunctionSpaceImp, class InducingFunctionImp>
class EllipticFiniteElement
{
public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef InducingFunctionImp InducingFunctionType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef Dune::Functionals::Common::LocalMatrix<RangeFieldType> LocalMatrixType;

  typedef Dune::Functionals::Common::LocalVector<RangeFieldType> LocalVectorType;

  typedef Dune::Functionals::Functional::L2<DiscreteFunctionSpaceImp, InducingFunctionImp> FunctionalType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  typedef typename DiscreteFunctionSpaceType::IteratorType EntityIteratorType;

  typedef typename EntityIteratorType::Entity EntityType;

  typedef typename EntityType::Geometry EntityGeometryType;

  typedef CachingQuadrature<GridPartType, 0> EntityQuadratureType;

  EllipticFiniteElement(DiscreteFunctionSpaceType& discreteFunctionSpace, InducingFunctionType& inducingFunction)
    : discreteFunctionSpace_(discreteFunctionSpace)
    , inducingFunction_(inducingFunction)
  {
  }

  ~EllipticFiniteElement()
  {
  }

  DiscreteFunctionSpaceType& space() const
  {
    return discreteFunctionSpace_;
  }

  template <class DiscreteFunctionImp1>
  const FunctionalType operator()(const DiscreteFunctionImp1& function1) const
  {
    FunctionalType func;
    return func;
  }

  template <class DiscreteFunctionImp1, class DiscreteFunctionImp2>
  const RangeFieldType operator()(const DiscreteFunctionImp1& function1, const DiscreteFunctionImp2& function2) const
  {
    return 0.0;
  }

  //  template< class DiscreteFucntionImp1, class DiscreteFunctionImp2 >
  //  const RangeFieldType operator()(  const DiscreteFucntionImp1& function1,
  //                                    const DiscreteFunctionImp2& function2 ) const
  //  {
  //    RangeFieldType ret = 0.0;

  //    // some types we will need
  //    typedef DiscreteFunctionImp1
  //      DiscreteFunctionType1;
  //    typedef DiscreteFunctionImp2
  //      DiscreteFunctionType2;

  //    typedef typename FunctionSpaceType::JacobianRangeType
  //      JacobianRangeType;

  //    // do gridwalk
  //    const EntityIteratorType BehindLastEntity = discreteFunctionSpace_.end();
  //    for ( EntityIteratorType entityIterator = discreteFunctionSpace_.begin();
  //          entityIterator != BehindLastEntity;
  //          ++entityIterator )
  //    {
  //      // some types we will need
  //      typedef typename EntityType::Geometry
  //        EntityGeometryType;
  //      typedef typename BaseFunctionSetType::RangeType
  //        RangeType;

  //      // entity and geometry
  //      const EntityType& entity = *entityIterator;
  //      const EntityGeometryType& entityGeometry = entity.geometry();

  //      // quadrature
  //      typedef CachingQuadrature< GridPartType, 0 >
  //        EntityQuadratureType;
  //      typedef typename EntityQuadratureType::CoordinateType
  //        EntityCoordinateType;
  //      const unsigned quadratureOrder = (2 * discreteFunctionSpace_.order() + 1);
  //      const EntityQuadratureType entityQuadrature( entity, quadratureOrder );
  //      const unsigned numberOfQuadraturePoints = entityQuadrature.nop();

  //      // do walk over quadrature points
  //      for(  unsigned int quadraturePoint = 0;
  //            quadraturePoint < numberOfQuadraturePoints;
  //            ++quadraturePoint )
  //      {
  //        // coordinates
  //        const EntityCoordinateType xReferenceElement = entityQuadrature.point( quadraturePoint );
  ////        const EntityCoordinateType xWorld = entityGeometry.global( xReferenceElement );
  //        const typename FunctionImp::DomainType xWorld = entityGeometry.global( xReferenceElement );

  //        // integration factors
  //        const double integrationFactor = entityGeometry.integrationElement( xReferenceElement );
  //        const double quadratureWeight = entityQuadrature.weight( quadraturePoint );

  //        // evaluate gradients of the analytical functions
  //        typename FunctionSpaceType::RangeType gradient1( 0.0 );
  //        JacobianRangeType gradient2( 0.0 );
  //        function1.evaluate( xWorld, gradient1 );
  ////        function2.jacobian( xWorld, gradient2 );

  //        // evaluate the inducing function
  //        RangeType functionValue = 0.0;
  //        inducingFunction_.evaluate( xWorld, functionValue );

  //        // compute integral
  ////        ret += integrationFactor * quadratureWeight * functionValue * gradient1 * gradient1;

  //      } // done walk over quadrature points


  //    } // done gridwalk

  //    return ret;
  //  }


  LocalMatrixType applyLocal(const EntityType& entity) const
  {

    // basefunctionset
    const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace_.baseFunctionSet(entity);
    const unsigned numberOfLocalDoFs          = baseFunctionSet.numBaseFunctions();

    // init return matrix
    LocalMatrixType ret(numberOfLocalDoFs, numberOfLocalDoFs);

    // geometry
    const EntityGeometryType& entityGeometry = entity.geometry();

    // quadrature
    typedef typename EntityQuadratureType::CoordinateType EntityCoordinateType;
    const unsigned quadratureOrder = discreteFunctionSpace_.order();
    const EntityQuadratureType entityQuadrature(entity, quadratureOrder);
    const unsigned numberOfQuadraturePoints = entityQuadrature.nop();

    // do loop over all local DoFs (i)
    for (unsigned int i = 0; i < numberOfLocalDoFs; ++i) {

      // do loop over all local DoFs (i)
      for (unsigned int j = 0; j < numberOfLocalDoFs; ++j) {

        // value of the operator, applied to the local basefunctions
        RangeFieldType operator_i_j = 0.0;

        // do walk over quadrature points
        for (unsigned int quadraturePoint = 0; quadraturePoint < numberOfQuadraturePoints; ++quadraturePoint) {
          // coordinates
          const EntityCoordinateType xReferenceElement = entityQuadrature.point(quadraturePoint);
          const EntityCoordinateType xWorld            = entityGeometry.global(xReferenceElement);

          // integration factors
          const double integrationFactor = entityGeometry.integrationElement(xReferenceElement);
          const double quadratureWeight  = entityQuadrature.weight(quadraturePoint);

          // jacobian inverse transposed
          const typename EntityGeometryType::Jacobian& jacobianInverseTransposed =
              entityGeometry.jacobianInverseTransposed(xReferenceElement);

          // evaluate function and basefunction
          RangeType functionValue = 0.0;
          inducingFunction_.evaluate(xWorld, functionValue);

          JacobianRangeType gradientBaseFunction_i_untransposed(0.0);
          baseFunctionSet.jacobian(i, xReferenceElement, gradientBaseFunction_i_untransposed);
          JacobianRangeType gradientBaseFunction_i(0.0);
          jacobianInverseTransposed.mv(gradientBaseFunction_i_untransposed[0], gradientBaseFunction_i[0]);

          JacobianRangeType gradientBaseFunction_j_untransposed(0.0);
          baseFunctionSet.jacobian(i, xReferenceElement, gradientBaseFunction_j_untransposed);
          JacobianRangeType gradientBaseFunction_j(0.0);
          jacobianInverseTransposed.mv(gradientBaseFunction_j_untransposed[0], gradientBaseFunction_j[0]);

          const double product = gradientBaseFunction_j[0] * gradientBaseFunction_j[0];

          // compute integral
          operator_i_j += integrationFactor * quadratureWeight * functionValue * product;

        } // done walk over quadrature points

        // set local matrix
        ret.set(i, j, operator_i_j);

      } // done loop over all local DoFs (i)

    } // done loop over all local DoFs (i)
  }

private:
  DiscreteFunctionSpaceType& discreteFunctionSpace_;
  InducingFunctionType& inducingFunction_;
};

} // end of namespace Operator

} // end of namespace Functionals

} // end of namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_OPERATOR_ELLIPTICFINITEELEMENT_HH
