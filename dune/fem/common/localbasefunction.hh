#ifndef DUNE_FEM_FUNCTIONALS_COMMON_LOCALBASEFUNCTIONSET_HH
#define DUNE_FEM_FUNCTIONALS_COMMON_LOCALBASEFUNCTIONSET_HH

// dune-fem includes
#include <dune/fem/function/localfunction/localfunction.hh>


namespace Dune
{

namespace Functionals
{

//! Contains several common classes.
namespace Common
{

/**
 * @brief Provides the selection of a local basis function out of a 
 * basis function set depending on the entity.
 *
 * Let @f$\tau_h@f$ be the triangulation of our domain @f$\Omega@f$.
 * Then for each entity @f$e@f$ in our grid there is set of basis functions 
 * @f$\varphi_e:=\{\varphi_{e,1},\ldots,\varphi_{e,n_e}\}@f$, where @f$n_e@f$ is the number of
 * local dof numbers depending on the entity @f$e@f$.
 *
 * Together with the provide() method you can select a single basis function 
 * @f$\varphi_{e',i}@f$ for an entity @f$e'@f$ and a local dof number @f$i@f$
 * from this basis function set.
 *
 * @tparam DiscreteFunctionSpaceImp The discrete function space type.
 */
template< class DiscreteFunctionSpaceImp >
class LocalBaseFunctionProvider
{
public:

  /**
   * @brief Represents an local basis function and all necessary methods.
   */
  class LocalBaseFunction
  {

  public:

    //! Type of the discret function space type.
    typedef DiscreteFunctionSpaceImp
      DiscreteFunctionSpaceType;

    //! Intrinsic type used for values in the domain field (usually a double)
    typedef typename DiscreteFunctionSpaceType::DomainFieldType
      DomainFieldType;

    //! Intrinsic type used for values in the range field (usually a double)
    typedef typename DiscreteFunctionSpaceType::RangeFieldType
      RangeFieldType;

    //! Type of domain vector (using type of domain field) has a Dune::FieldVector type interface.
    typedef typename DiscreteFunctionSpaceType::DomainType
      DomainType;

    //! Type of range vector (using type of range field) has a Dune::FieldVector type interface.
    typedef typename DiscreteFunctionSpaceType::RangeType
      RangeType;

    //! Intrinsic type used for the jacobian values has a Dune::FieldMatrix type interface.
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType
      JacobianRangeType;

    //! Intrinsic type used for the hessian values has a Dune::FieldMatrix type interface.
    typedef typename DiscreteFunctionSpaceType::HessianRangeType
      HessianRangeType;

    //! Type of the base function set.
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
      BaseFunctionSetType;

    //! Type of the grid part.
    typedef typename DiscreteFunctionSpaceType::GridPartType
      GridPartType;

    //! Type of the iterator iterating over entities.
    typedef typename DiscreteFunctionSpaceType::IteratorType
      EntityIteratorType;

    //! Type of entities.
    typedef typename EntityIteratorType::Entity
      EntityType;

    //! Type of geometry of the entities.
    typedef typename EntityType::Geometry
      EntityGeometryType;

    //! Dimension of the domain.
    static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

    //! Dimension of the range.
    static const int dimRange = DiscreteFunctionSpaceType::dimRange;

    /**
     * @brief Constructor, selecting information about the basis function @f$\varphi_{e',i}@f$
     * which can be extracted from the basis function set @f$\varphi_{e'}@f$.
     *
     * @param entity The entity @f$e'@f$ the basefunction set @f$\varphi_{e'}@f$ belongs to.
     * @param baseFunctionSet The basisfunction set @f$\varphi_{e'}@f$.
     * @param localDoFNumber The local dof number.
     */
    LocalBaseFunction( const EntityType& entity,
                       const BaseFunctionSetType& baseFunctionSet,
                       const int localDoFNumber )
      : entity_( entity ),
        baseFunctionSet_( baseFunctionSet ),
        localDoFNumber_( localDoFNumber )
    {
    }

    /**
     * @brief Deconstructor.
     */
    ~LocalBaseFunction()
    {
    }

    /**
     * @brief Returns the polynomial order of the discrete function space.
     *
     * @return The polynomial order of the discrete function space.
     */
    int order() const
    {
      return DiscreteFunctionSpaceType::polynomialOrder;
    }

    /**
     * @brief Evaluates a basis function @f$\varphi_{e',i}@f$ on a given
     * local point @f$x@f$.
     *
     * In short: Computes @f$\varphi_{e',i}(x)@f$.
     *
     * @param x The point @f$x@f$, where the basis function 
     * @f$\varphi_{e',i}@f$ should be evaluated.
     * @param ret The result of the evaluation, i.e. @f$\varphi_{e',i}(x)@f$.
     */
    template< class PointType >
    void evaluate( const PointType& x, RangeType& ret) const
    {
      baseFunctionSet_.evaluate( localDoFNumber_, x, ret );
    }

    /**
     * @brief Evaluates the jacobian of the basefunction set on a given
     * local point @f$x@f$.
     *
     * In short: Computes @f$\nabla\varphi_{e',i}(x)@f$.
     *
     * @param x The point @f$x@f$, where the gradient of the basis 
     * function @f$\nabla\varphi_{e',i}@f$ should be evaluated.
     * @param ret The result of the evaluation of the jacobian of the basis function,
     * i.e. @f$\nabla\varphi_{e',i}(x)@f$.
     */
    template< class PointType >
    void jacobian( const PointType& x, JacobianRangeType& ret ) const
    {
      // some types we will need
      typedef typename EntityGeometryType::Jacobian
        JacobianInverseTransposedType;

      typedef typename JacobianRangeType::row_type
        JacobianRowType;

      // geometry and jacobian inverse transposed
      const EntityGeometryType& entityGeometry = entity_.geometry();
      const JacobianInverseTransposedType& jacobianInverseTransposed = entityGeometry.jacobianInverseTransposed( x );

      // get untransposed jacobian
      JacobianRangeType jacobianUntransposed( 0.0 );
      baseFunctionSet_.jacobian( localDoFNumber_, x, jacobianUntransposed );

      // do for each dim of Range
      for( unsigned int row = 0; row < ret.N(); ++row )
      {
        // transpose
        JacobianRowType jacobian( 0.0 );
        jacobianInverseTransposed.mv( jacobianUntransposed[0], jacobian );

        // return
        ret[row] = jacobian;
      }
    }

    /**
     * @brief Returns a reference to the entity @f$e'@f$ the basis function @f$\varphi_{e',i}@f$ is living on.
     *
     * @return A reference to the entity @f$e'@f$ the basis function @f$\varphi_{e',i}@f$ is living on.
     */
    const EntityType& entity() const
    {
      return entity_;
    }

  private:

    const EntityType& entity_;
    const BaseFunctionSetType baseFunctionSet_;
    const int localDoFNumber_;

  }; // end class LocalBaseFunction

  //! Type of the local base functions @f$\varphi_{e',i}@f$.
  typedef LocalBaseFunction
    LocalBaseFunctionType;

  //! Type of the discrete function space of the local basis function.
  typedef typename LocalBaseFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  //! Type of the entity @f$e@f$ the basefunction is living on.
  typedef typename LocalBaseFunctionType::EntityType
    EntityType;

  /**
   * @brief Constructor.
   *
   * @param discreteFunctionSpace The discrete function space.
   */
  LocalBaseFunctionProvider( const DiscreteFunctionSpaceType& discreteFunctionSpace )
    : discreteFunctionSpace_( discreteFunctionSpace )
  {
  }

  /**
   * @brief Deconstructor.
   */
  ~LocalBaseFunctionProvider()
  {
  }

  /**
   * @brief Returns a local base function for a given entity and a local dof number.
   *
   * @param entity The entity the base function is living on.
   * @param localDofNumber The number of the base function we want to evaluate.
   */
  const LocalBaseFunctionType provide( const EntityType& entity, const int localDoFNumber ) const
  {
    return LocalBaseFunctionType( entity, discreteFunctionSpace_.baseFunctionSet( entity ), localDoFNumber );
  }

private:

  const DiscreteFunctionSpaceType& discreteFunctionSpace_;

}; // end class LocalBaseFunctionProvider

} // end namespace Common

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_COMMON_LOCALBASEFUNCTIONSET_HH
