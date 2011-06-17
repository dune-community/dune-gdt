#ifndef DUNE_FUNCTIONALS_COMMON_LOCALBASEFUNCTION_HH
#define DUNE_FUNCTIONALS_COMMON_LOCALBASEFUNCTION_HH


namespace Dune
{

namespace Functionals
{

//! Contains several common classes.
namespace Common
{

template< class DiscreteFunctionSpaceImp >
class LocalBaseFunctionSet
{
private:

  class LocalBaseFunction
  {

  public:

    //! Type of the discret function space type.
    typedef DiscreteFunctionSpaceImp
      DiscreteFunctionSpaceType;

    typedef LocalBaseFunctionSet< DiscreteFunctionSpaceType >
      LocalBaseFunctionSetType;

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

    //! Type of entities.
    typedef typename DiscreteFunctionSpaceType::EntityType
      EntityType;

    /**
     * @brief Constructor, selecting information about the basis function @f$\varphi_{e',i}@f$
     * which can be extracted from the basis function set @f$\varphi_{e'}@f$.
     *
     * @param entity The entity @f$e'@f$ the basefunction set @f$\varphi_{e'}@f$ belongs to.
     * @param baseFunctionSet The basisfunction set @f$\varphi_{e'}@f$.
     * @param localDoFNumber The local dof number.
     */
    LocalBaseFunction( const LocalBaseFunctionSetType& localBaseFunctionSet,
                       const int localDoFNumber )
      : localBaseFunctionSet_( localBaseFunctionSet ),
        entity_( localBaseFunctionSet.entity() ),
        localDoFNumber_( localDoFNumber )
    {
    }

    /**
     * @brief Deconstructor.
     */
    ~LocalBaseFunction()
    {
    }

    const LocalBaseFunctionSetType& baseFunctionSet() const
    {
      return localBaseFunctionSet_;
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

    const int number() const
    {
      return localDoFNumber_;
    }

    /**
     * @brief Returns the polynomial order of the discrete function space.
     *
     * @return The polynomial order of the discrete function space.
     */
    int order() const
    {
      return localBaseFunctionSet_.order();
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
    void evaluate( const DomainType& x, RangeType& ret) const
    {
      localBaseFunctionSet_.evaluate( localDoFNumber_, x, ret );
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
    void jacobian( const DomainType& x, JacobianRangeType& ret ) const
    {
      localBaseFunctionSet_.jacobian( localDoFNumber_, x, ret );
    }

  private:

    const LocalBaseFunctionSetType& localBaseFunctionSet_;
    const EntityType& entity_;
    const int localDoFNumber_;

  }; // end class LocalBaseFunction

public:

  typedef LocalBaseFunction
    LocalBaseFunctionType;

  typedef typename LocalBaseFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  typedef typename LocalBaseFunctionType::EntityType
    EntityType;

  typedef typename LocalBaseFunctionType::DomainFieldType
    DomainFieldType;

  typedef typename LocalBaseFunctionType::DomainType
    DomainType;

  typedef typename LocalBaseFunctionType::RangeFieldType
    RangeFieldType;

  typedef typename LocalBaseFunctionType::RangeType
    RangeType;

  typedef typename LocalBaseFunctionType::JacobianRangeType
    JacobianRangeType;

  typedef typename LocalBaseFunctionType::HessianRangeType
    HessianRangeType;

private:

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    HostBaseFunctionSetType;

public:

  LocalBaseFunctionSet( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    : space_( space ),
      entity_( entity ),
      hostBaseFunctionSet_( space.baseFunctionSet( entity ) )
  {
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  const LocalBaseFunctionType baseFunction( const int i ) const
  {
    return LocalBaseFunctionType( *this, i );
  }

  const int order() const
  {
    return space_.order();
  }

  const int numBaseFunctions() const
  {
    return hostBaseFunctionSet_.numBaseFunctions();
  }

  void evaluate( const int i, const DomainType& x, RangeType& ret ) const
  {
    hostBaseFunctionSet_.evaluate( i, x, ret );
  }

  /**
    \brief      evaluates the jacobian of the ith ocal basefunction
    \attention  the evalaution is already multiplied by entityGeometry.jacobianInverseTransposed( x )
    **/
  void jacobian( const int i, const DomainType& x, JacobianRangeType& ret ) const
  {
    // some types we will need
    typedef typename EntityType::Geometry
      EntityGeometryType;

    typedef typename EntityGeometryType::Jacobian
      JacobianInverseTransposedType;

    typedef typename JacobianRangeType::row_type
      JacobianRowType;

    // geometry and jacobian inverse transposed
    const EntityGeometryType& entityGeometry = entity_.geometry();
    const JacobianInverseTransposedType& jacobianInverseTransposed = entityGeometry.jacobianInverseTransposed( x );

    // get untransposed jacobian
    JacobianRangeType jacobianUntransposed( 0.0 );
    hostBaseFunctionSet_.jacobian( i, x, jacobianUntransposed );

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

private:

  const DiscreteFunctionSpaceType& space_;
  const EntityType& entity_;
  const HostBaseFunctionSetType hostBaseFunctionSet_;

}; // end class LocalBaseFunctionSet

} // end namespace Common

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_COMMON_LOCALBASEFUNCTION_HH
