#ifndef DUNE_FUNCTIONALS_EVALUATION_UNARY_PRODUCT_HH
#define DUNE_FUNCTIONALS_EVALUATION_UNARY_PRODUCT_HH

namespace Dune
{

namespace Functionals
{

namespace Evaluation
{

namespace Unary
{

/**
  \brief  This represents the operation \f$fv\f$.

          \f$f\f$ is a given right hand side (in this case 1) and \f$v\f$ may be a local function, i.e. a
          testfunction.
  \tparam FunctionSpaceImp
          Type of the function space, where \f$f\f$ and \f$v\f$ live in.
  **/
template< class FunctionSpaceImp >
class Product
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef Dune::FemTools::Function::Runtime< FunctionSpaceType >
    InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  Product( const std::string expression = "[1.0;1.0;1.0]" )
    : inducingFunction_( expression ),
      order_( 1 )
  {
  }

  //! constructor, takes the inducing functions expression as a runtime parameter
  Product( const std::string expression = "[1.0;1.0;1.0]", const int order = 1 )
    : inducingFunction_( expression ),
      order_( 1 )
  {
    if( order < 0 )
      order_ = 0;
    else
      order_ = order;
  }

  //! copy constructor
  Product( const Product& other )
    : inducingFunction_( other.inducingFunction() ),
      order_( other.order() )
  {
  }

  //! returns the inducing function
  const InducingFunctionType& inducingFunction() const
  {
    return inducingFunction_;
  }

  unsigned int order() const
  {
    return order_;
  }

  /**
    \brief      Evaluates \f$f(x)v(x)\f$ for a given local point \f$x\f$.

    \tparam     LocalTestFunctionType
                Type of the local function \f$v\f$, i.e. Dune::LocalFunction.
    \tparam     LocalPointType
                Type of the local point \f$x\f$, i.e. Dune::FieldVector.
    \param[in]  localTestFunction
                The local function \f$v\f$.
    \param[in]  localPoint
                The local point \f$x\f$. This point is local in the sense, that this is a point on a reference
                element.
    \return     \f$f(x)v(x)\f$
    **/
  template< class LocalTestBaseFunctionSetType, class LocalVectorType >
  void evaluate(  const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                  const DomainType& localPoint,
                  LocalVectorType& ret ) const
  {
    // get global point
    const DomainType globalPoint = localTestBaseFunctionSet.entity().geometry().global( localPoint );

    // evaluate inducing function
    RangeType functionValue( 0.0 );
    inducingFunction_.evaluate( globalPoint, functionValue );

    // evaluate set of local functions
    const unsigned int size = localTestBaseFunctionSet.size();
    std::vector< RangeType > valuesLocalBaseFunctionSet( size, RangeType( 0.0 ) );
    localTestBaseFunctionSet.evaluate( localPoint, valuesLocalBaseFunctionSet );

    // do loop over all basis functions
    assert( ret.size() == size );
    for( unsigned int i = 0; i < size; ++i )
    {
      ret[i] = functionValue * valuesLocalBaseFunctionSet[i];
    }
  }

private:

  const InducingFunctionType inducingFunction_;
  unsigned int order_;
}; // end class Product

} // end namespace Unary

} // end namespace Evaluation

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_EVALUATION_UNARY_PRODUCT_HH
