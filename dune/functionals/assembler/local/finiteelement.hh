#ifndef DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_FINITEELEMENT_HH
#define DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_FINITEELEMENT_HH

// dune-functionals includes
#include <dune/functionals/common/localmatrix.hh>

namespace Dune
{

namespace Functionals
{

namespace Assembler
{

namespace Local
{

template< class LocalOperatorImp >
class ContinuousFiniteElement
{
public:

  typedef LocalOperatorImp
    LocalOperatorType;

  typedef typename LocalOperatorType::RangeFieldType
    RangeFieldType;

  typedef Dune::Functionals::Common::LocalMatrix< RangeFieldType >
    LocalMatrixType;

  //! constructor
  ContinuousFiniteElement( const LocalOperatorType& localOperator )
    : localOperator_( localOperator )
  {
  }

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

  template< class LocalAnsatzBaseFunctionSetType, class LocalTestBaseFunctionSetType >
  void assembleLocal( const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                      const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                      LocalMatrixType& localMatrix ) const
  {
    localOperator_.applyLocal( localAnsatzBaseFunctionSet, localTestBaseFunctionSet, localMatrix );
  }

private:

  const LocalOperatorType& localOperator_;

}; // end class ContinuousFiniteElement

} // end namespace Local

} // end namespace Assembler

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_FINITEELEMENT_HH
