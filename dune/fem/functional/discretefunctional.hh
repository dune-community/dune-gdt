#ifndef DUNE_FEM_DISCRETEFUNCTIONAL_HH
#define DUNE_FEM_DISCRETEFUNCTIONAL_HH

namespace Dune
{

  namespace Fem
  {

    template< class DiscreteFunction >
    class DiscreteFunctional
    {
      typedef DiscreteFunctional< DiscreteFunction > ThisType;

    public:
      typedef DiscreteFunction DiscreteFunctionType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      DiscreteFunctional ( const std::string &name, const DiscreteFunctionSpace &dfSpace );

      FieldType operator () ( const DiscreteFunction &u ) const;

      const std::string &name () const;
      const DiscreteFunctionSpace &space () const;

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DISCRETEFUNCTIONAL_HH
