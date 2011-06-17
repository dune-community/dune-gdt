/**
  \file   integration.hh
  **/

#ifndef DUNE_FEM_FUNCTIONALS_DISCRETEOPERATOR_LOCAL_INTEGRATION_HH
#define DUNE_FEM_FUNCTIONALS_DISCRETEOPERATOR_LOCAL_INTEGRATION_HH

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-functionals includes
#include <dune/functionals/common/localmatrix.hh>
#include <dune/functionals/common/localvector.hh>

namespace Dune
{

//! dune-fem-functionals
namespace Functionals
{

//! Contains several operators.
namespace DiscreteOperator
{

namespace Local
{

/**
  \todo   correct doc
  \brief  Linear operator.

          This class represents a linear operator \f$A: V_{h} \rightarrow W_{h}^{-1}\f$ which arises from a weak
          formulation, where \f$V_{h}\f$ and \f$W_{h}\f$ are discrete function spaces, i.e. Lagrange finite element
          spaces. These are usually called ansatz and test space. A linear operator induces a bilinear
          form
          \f[
            A: V_h \times W_h \rightarrow \mathbb{R} \text{,}
          \f]
          \f[
            u,v \mapsto A(u)[v] \text{,}
          \f]
          where \f$A(u): W_h \rightarrow \mathbb{R}\f$ itself is a functional for each \f$u \in V_h\f$. Because of the
          linearity of \f$A\f$ we can decompose the the application of \f$A(u)[\psi]\f$ to an argument
          \f$u \in V_{h}\f$, testet by a basefunction \f$\psi \in W_h\f$ into
          local applications of \f$A\f$ to the local basefunctions of the ansatz and test space:
          \f[
            A(u)[\psi] = \sum_{E \in \mathcal{T}_{h}} \sum_{i \in I_{E}} \sum_{j \in J_{E}} u_i A( \varphi_i )[\psi_j]\text{,}
          \f]
          where \f$E\f$ are the codim 0 entities of a triangulation \f$\mathcal{T}_{h}\f$, \f$I_E\f$ and \f$J_E\f$ are
          index sets of local DoFs, \f$u_{i}\f$ are the corresponding coefficients, \f$\varphi_i\f$ are the
          local basefunctions of the ansatz space \f$V_{h}\f$ and \f$\psi_j\f$ are the local basefunctions of the
          test space \f$W_h\f$. A linear operator thus only hast to provide this local application
          \f$A( \varphi_i )[\psi_j]\f$. This local application (and thus the linear operator) is induced by a local
          operation \f$a(\varphi,\psi)\f$ which operates on local functions. In the finite element case for example,
          where the operator is given by
          \f[
            A(u)[v] = \int_{\Omega} a \nabla \varphi \nabla \psi \text{dx} \text{,}
          \f]
          for some function \f$a\f$, the corresponding local operation is given by
          \f[
            a(\varphi,\psi) := \int_{E} a(x) \nabla \varphi(x) \nabla \psi(x) \text{dx} \text{.}
          \f]
  \tparam LocalOperationImp
          Type of the local operation \f$a\f$, i.e. Dune::Functionals::LocalOperation::Interface. This local operation
          has to implement the method operate( localTestFunction, localAnsatzFunction ).
  \tparam DiscreteAnsatzFunctionSpaceImp
          Type of the discrete function space \f$V_h\f$, where the ansatz functions live in, i.e.
          Dune::DiscreteFunctionSpaceInterface.
  \tparam DiscreteTestFunctionSpaceImp
          Type of the discrete function space \f$W_h\f$, where the test functions live in, i.e.
          Dune::DiscreteFunctionSpaceInterface. If the test and the ansatz space are identicall, this template
          argument may be omitted.
  **/
template< class LocalEvaluationImp >
class Codim0Integration
{
public:

  typedef LocalEvaluationImp
    LocalEvaluationType;

  typedef LocalEvaluationType
    FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef Dune::Functionals::Common::LocalMatrix< RangeFieldType >
    LocalMatrixType;

  /**
    \brief      Constructor storing the local operation and the ansatz and test space.

                Use this constructor, if ansatz and test space are not the same.
    \param[in]  localOperation
                The local operation \f$a\f$, which induces the operator.
    \param[in]  ansatzSpace
                The space of ansatz functions \f$V_h\f$.
    \param[in]  testSpace
                The space of test functions \f$W_h\f$.
    **/
  Codim0Integration(  const LocalEvaluationType& localEvaluation )
    : localEvaluation_( localEvaluation )
  {
  }

  /**
    \brief  Returns the local operation \f$a\f$.
    \return \f$a\f$.
    **/
  const LocalEvaluationType localEvaluation() const
  {
    return localEvaluation_;
  }

  /**
    \brief      Local application of the operator.

                This method represents the application of the operator to all local basefunctions of the ansatz and
                test space on a given entity:
                \f[
                  \{A( \varphi_i )[\psi_j]\}_{i \in I_E, j \in J_E} \text{.}
                \f]
    \param[in]  entity
                The entity, on wich the operator is being applied on.
    \return     The matrix \f$\{A( \varphi_i )[\psi_j]\}_{i \in I_E, j \in J_E}\f$.
    **/
  template< class LocalAnsatzBaseFunctionSetType, class LocalTestBaseFunctionSetType >
  void applyLocal( const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                   const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                   LocalMatrixType& localMatrix ) const
  {
    // some types
    typedef typename LocalAnsatzBaseFunctionSetType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType
      GridPartType;

    typedef Dune::CachingQuadrature< GridPartType, 0 >
      VolumeQuadratureType;

    typedef typename LocalAnsatzBaseFunctionSetType::LocalBaseFunctionType
      LocalAnsatzBaseFunctionType;

    typedef typename LocalTestBaseFunctionSetType::LocalBaseFunctionType
      LocalTestBaseFunctionType;

    // some stuff
    const unsigned numberOfLocalAnsatzDoFs = localAnsatzBaseFunctionSet.numBaseFunctions();
    const unsigned numberOfLocalTestDoFs = localTestBaseFunctionSet.numBaseFunctions();
    const unsigned int quadratureOrder = 1 + localAnsatzBaseFunctionSet.order() + localTestBaseFunctionSet.order();
    const VolumeQuadratureType volumeQuadrature( localAnsatzBaseFunctionSet.entity(), quadratureOrder );
    const unsigned int numberOfQuadraturePoints = volumeQuadrature.nop();

    // do loop over all local ansatz DoFs
    for( unsigned int i = 0; i < numberOfLocalAnsatzDoFs; ++i )
    {
      const LocalAnsatzBaseFunctionType localAnsatzBaseFunction_i = localAnsatzBaseFunctionSet.baseFunction( i );

      // do loop over all local test DoFs
      for( unsigned int j = 0; j < numberOfLocalTestDoFs; ++j )
      {
        const LocalTestBaseFunctionType localTestBaseFunction_j = localTestBaseFunctionSet.baseFunction( j );

        // do loop over all quadrature points
        RangeFieldType operator_i_j( 0.0 );
        for( unsigned int q = 0; q < numberOfQuadraturePoints; ++q )
        {
          // local coordinate
          const DomainType x = volumeQuadrature.point( q );

          // integration factors
          const double integrationFactor = localAnsatzBaseFunctionSet.entity().geometry().integrationElement( x );
          const double quadratureWeight = volumeQuadrature.weight( q );

          // evaluate the local operation
          const RangeFieldType localOperationEvalauted = localEvaluation_.evaluate( localAnsatzBaseFunction_i, localTestBaseFunction_j, x );

          // compute integral
          operator_i_j += integrationFactor * quadratureWeight * localOperationEvalauted;
        } // done loop over all quadrature points

        // set local matrix (the = is important, since we dont assume a clean matrix)
        localMatrix[i][j] = operator_i_j;

      } // done loop over all local test DoFs

    } // done loop over all local ansatz DoFs

  } // end method applyLocal

private:

  const LocalEvaluationType& localEvaluation_;

}; // end class Codim0Integration

} // end namespace Local

} // end namespace DiscreteOperator

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_DISCRETEOPERATOR_LOCAL_INTEGRATION_HH
