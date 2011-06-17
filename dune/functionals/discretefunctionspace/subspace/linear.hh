#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH

// dune-functionals includes
#include <dune/functionals/constraints/dirichlet.hh>

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunctionSpace
{

namespace Subspace
{

template< class SuperSpaceImp >
class DirichletZero
{
public:

  typedef SuperSpaceImp
    SuperSpaceType;

  typedef Dune::Functionals::Constraints::DirichletZero< SuperSpaceType >
    DirichletZeroConstraintsType;

  DirichletZero( const SuperSpaceType& superSpace )
    : superSpace_( superSpace ),
      constraints_( superSpace_ )
  {
  }

//  /**
//   * @brief Copy constructor.
//   *
//   * @param lin A reference to an existing linear subspace @f$V'_{h,C}@f$.
//   */
//  Linear( const Linear& lin )
//    : DiscreteFunctionSpaceType( lin.space().gridPart() ),
//      space_( lin.space() ),
//      constraints_( lin.constraints() )
//  {
//  }

//  /**
//   * @brief Returns the constraints.
//   *
//   * @return A reference to the constraints @f$C@f$.
//   */
//  ConstraintsType& constraints() const
//  {
//    return constraints_;
//  }

//  /**
//   * @brief Returns the discrete function space.
//   *
//   * @return A reference to the discrete function space @f$V_h@f$.
//   */
//  DiscreteFunctionSpaceType& space() const
//  {
//    return space_;
//  }

private:
  const SuperSpaceType& superSpace_;
  const DirichletZeroConstraintsType constraints_;

}; // end of class DirichletZero

} // end of namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH */
