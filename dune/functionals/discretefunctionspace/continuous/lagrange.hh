#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>

// dune-pdelab includes
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

// dune-functionals includes
#include <dune/functionals/common/localbasefunctionset.hh>

namespace Dune {

namespace Functionals {

namespace DiscreteFunctionSpace {

namespace Continuous {

template <class FunctionSpaceImp, class GridViewImp, int polOrder>
class Lagrange
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef GridViewImp GridViewType;

  typedef Lagrange<FunctionSpaceType, GridViewType, polOrder> ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  static const unsigned int dimDomain = FunctionSpaceType::dimDomain;

  static const unsigned int dimRange = FunctionSpaceType::dimRange;

private:
  typedef Dune::PDELab::P1LocalFiniteElementMap<DomainFieldType, RangeFieldType, dimDomain> LocalFiniteElementMapType;

public:
  //  typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  //    HostSpaceType;

  typedef Dune::PDELab::GridFunctionSpace<GridViewType, LocalFiniteElementMapType> HostSpaceType;

  typedef Dune::Functionals::Common::LocalBaseFunctionSet<ThisType> LocalBaseFunctionSetType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;

  //  /**
  //    \defgroup dune-fem related
  //    \{
  //    **/
  //  typedef typename HostSpaceType::BaseFunctionSetType
  //    BaseFunctionSetType;

  //  typedef typename HostSpaceType::IteratorType
  //    IteratorType;
  //  /**
  //    \}
  //    **/


  Lagrange(GridViewType& gridView)
    : gridView_(gridView)
    , localFiniteElementMap_()
    , hostSpace_(gridView_, localFiniteElementMap_) /*,
     numMaxLocalDoFs_( -1 )*/
  {
    //    // in the simple case, there should be the same number of dofs on each entity
    //    const IteratorType entityIterator = begin();
    //    const EntityType& entity = *entityIterator;
    //    numMaxLocalDoFs_ = hostSpace_.baseFunctionSet( entity ).numBaseFunctions();
  }

  //  const GridPartType& gridPart() const
  //  {
  //    return gridPart_;
  //  }

  //  const HostSpaceType& hostSpace() const
  //  {
  //    return hostSpace_;
  //  }

  //  const LocalBaseFunctionSetType localBaseFunctionSet( const EntityType& entity ) const
  //  {
  //    return LocalBaseFunctionSetType( *this, entity );
  //  }

  const unsigned int size() const
  {
    return hostSpace_.size();
  }

  //  const int numMaxLocalDoFs() const
  //  {
  //    return numMaxLocalDoFs_;
  //  }

  //  const int order() const
  //  {
  //    return hostSpace_.order();
  //  }

  //  /**
  //    \defgroup dune-fem related
  //    \{
  //    **/
  //  IteratorType begin() const
  //  {
  //    return hostSpace_.begin();
  //  }

  //  const IteratorType end() const
  //  {
  //    return hostSpace_.end();
  //  }

  //  const BaseFunctionSetType baseFunctionSet( const EntityType& entity ) const
  //  {
  //    return hostSpace_.baseFunctionSet( entity );
  //  }

  //  int mapToGlobal( const EntityType& entity, const int localDof) const
  //  {
  //    return hostSpace_.mapToGlobal( entity, localDof);
  //  }
  //  /**
  //    \}
  //    **/

private:
  const GridViewType& gridView_;
  const LocalFiniteElementMapType localFiniteElementMap_;
  const HostSpaceType hostSpace_;
  //  unsigned int numMaxLocalDoFs_;

}; // end class Lagrange

} // end namespace Continuous

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGE_HH
