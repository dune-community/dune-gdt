#ifndef DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH
#define DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH


namespace Dune
{

namespace Functionals
{

namespace Subspace
{

// \todo should prepare sparsity patterns and such things!
template< class DiscreteFunctionSpaceImp, class ConstraintsImp >
class Linear
 : public DiscreteFunctionSpaceImp
{
public:
  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;
  typedef ConstraintsImp
    ConstraintsType;

  typedef typename DiscreteFunctionSpaceType::Traits
    Traits;
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType
    FunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType
    GridType;
  typedef typename DiscreteFunctionSpaceType::IndexSetType
    IndexSetType;
  typedef typename DiscreteFunctionSpaceType::IteratorType
    IteratorType;
  typedef typename DiscreteFunctionSpaceType::EntityType
    EntityType;
  typedef typename DiscreteFunctionSpaceType::DofManagerType
    DofManagerType;
  typedef typename DiscreteFunctionSpaceType::CommunicationManagerType
    CommunicationManagerType;
  typedef typename DiscreteFunctionSpaceType::MapperType
    MapperType;
  typedef typename DiscreteFunctionSpaceType::BlockMapperType
    BlockMapperType;
    
private:
  typedef DiscreteFunctionSpaceType
    BaseType;
  typedef typename FunctionSpaceType::DomainType
    DomainType;
  typedef typename FunctionSpaceType::RangeType
    RangeType;
  typedef typename FunctionSpaceType::DomainFieldType
    DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;


public:
  enum {localBlockSize = DiscreteFunctionSpaceType::localBlockSize};

  Linear( DiscreteFunctionSpaceType& space,
          ConstraintsType& constraints )
    : DiscreteFunctionSpaceType( space.gridPart() ),
      space_( space ),
      constraints_( constraints )
  {
  }

  Linear( const Linear& lin )
    : DiscreteFunctionSpaceType( lin.space().gridPart() ),
      space_( lin.space() ),
      constraints_( lin.constraints() )
  {
  }
  
  ConstraintsType& constraints() const
  {
    return constraints_;
  }

  DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

private:
  DiscreteFunctionSpaceType& space_;
  ConstraintsType& constraints_;
}; // end of class Linear

template< class LinearSubspaceImp, class OffsetFunctionImp >
class Affine
  : public LinearSubspaceImp
{
public:
  typedef LinearSubspaceImp
    LinearSubspaceType;    
  typedef OffsetFunctionImp
    OffsetFunctionType;
  
  typedef typename LinearSubspaceType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename LinearSubspaceType::ConstraintsType
    ConstraintsType;

  typedef typename LinearSubspaceType::Traits
    Traits;

  typedef typename LinearSubspaceType::FunctionSpaceType
    FunctionSpaceType;
  typedef typename LinearSubspaceType::BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename LinearSubspaceType::GridPartType
    GridPartType;
  typedef typename LinearSubspaceType::GridType
    GridType;
  typedef typename LinearSubspaceType::IndexSetType
    IndexSetType;
  typedef typename LinearSubspaceType::IteratorType
    IteratorType;
  typedef typename GridPartType::IntersectionIteratorType
    IntersectionIteratorType;
  typedef typename LinearSubspaceType::EntityType
    EntityType;
  typedef typename LinearSubspaceType::DofManagerType
    DofManagerType;
  typedef typename LinearSubspaceType::CommunicationManagerType
    CommunicationManagerType;
  typedef typename LinearSubspaceType::MapperType
    MapperType;
  typedef typename LinearSubspaceType::BlockMapperType
    BlockMapperType;

private:
  typedef LinearSubspaceImp
    BaseType;
  typedef typename FunctionSpaceType::DomainType
    DomainType;
  typedef typename FunctionSpaceType::RangeType
    RangeType;
  typedef typename FunctionSpaceType::DomainFieldType
    DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;


public:
  enum {localBlockSize = LinearSubspaceType::localBlockSize};

  Affine( LinearSubspaceType& linear,
          OffsetFunctionType& offset )
    : LinearSubspaceType( linear ),
      linear_( linear ),
      offset_( offset )
  {
  }

  Affine( const Affine& aff )
    : LinearSubspaceType( aff.linearSpace() ),
      linear_( aff.linearSpace() ),
      offset_( aff.offset() )
  {
  }
  
  LinearSubspaceType& linearSpace() const
  {
    return linear_;
  }

  OffsetFunctionType& offset() const
  {
    return offset_;
  }

private:
  LinearSubspaceType& linear_;
  OffsetFunctionType& offset_;
}; // end of class Affine

} // end of namespace Subspace

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH */
