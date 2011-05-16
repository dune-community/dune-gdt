#ifndef DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_INTERFACE_HH
#define DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_INTERFACE_HH

namespace Dune {

namespace Functionals {

namespace DiscreteFunctional {

template <class DiscreteFunction>
class Interface
{
  typedef Interface<DiscreteFunction> ThisType;

public:
  typedef DiscreteFunction DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  Interface(const std::string& name, const DiscreteFunctionSpace& dfSpace);

  FieldType operator()(const DiscreteFunction& u) const;

  const std::string& name() const;
  const DiscreteFunctionSpace& space() const;
};

} // end namespace DiscreteFunctional

} // end namespace Functionals

} // end namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTIONALS_DISCRETEFUNCTIONAL_INTERFACE_HH
