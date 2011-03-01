#ifndef DUNE_FEM_FUNCTIONALS_DISCRETELINEARFUNCTIONAL_HH
#define DUNE_FEM_FUNCTIONALS_DISCRETELINEARFUNCTIONAL_HH

// dune fem-functionals includes
#include <dune/fem/dofvector/dofvector.hh>

// include this file after all other includes because some of them might undef the macros we want to use
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

namespace Functionals
{

/**
  * \brief      This class is the interface for discrete linear functionals.
  *
  *             This class is mainly here for documentational purpose.
  *
  * \attention  This class is under construction!
  *
  * \todo       Doc me, please!
  **/
template< class DiscreteLinearFunctionalTraitsImp >
class DiscreteLinearFunctionalInterface
{
public:

  //! Traits.
  typedef DiscreteLinearFunctionalTraitsImp
    Traits;

private:

  //! For crtp trick.
  typedef typename Traits::DiscreteLinearFunctionalType
    DiscreteLinearFunctionalType;

public:

  //! Type of the function, which induces the functional.
  typedef typename Traits::InducingFunctionType
    InducingFunctionType;

  //! Type of the inducing functions range field.
  typedef typename InducingFunctionType::RangeFieldType
    RangeFieldType;

  //! Type of the local DoF vector.
  typedef Dune::Functionals::LocalDoFVector< RangeFieldType >
    LocalDoFVectorType;

  //! Constructor (empty).
  DiscreteLinearFunctionalInterface()
  {
    std::cout << "DiscreteLinearFunctionalInterface::DiscreteLinearFunctionalInterface()" << std::endl;
  }

  //! Destructor (empty).
  ~DiscreteLinearFunctionalInterface()
  {
    std::cout << "DiscreteLinearFunctionalInterface::~DiscreteLinearFunctionalInterface()" << std::endl;
  }

  /**
    * \brief      This operator represents the application of the functional to a discrete function.
    *
    *             This operator calls the operator of the derived class.
    *
    * \todo       Doc me, please!
    **/
  template< class DiscreteFunctionType >
  const RangeFieldType operator()( const DiscreteFunctionType& discreteFunction ) const
  {
    std::cout << "DiscreteLinearFunctionalInterface::operator()" << std::endl;
    CHECK_INTERFACE_IMPLEMENTATION( asImp().operator()( discreteFunction ) );
    return asImp().operator()( discreteFunction );
  }

  /**
    * \brief      This function represents the application of the functional to a local basefunctionset.
    *
    *             This function calls the function of the derived class.
    *
    * \todo       Doc me, please!
    **/
  template< class DiscreteFunctionSpaceType, class EntityType, class BaseFunctionSetType >
  const LocalDoFVectorType applyLocal( const DiscreteFunctionSpaceType& discreteFunctionSpace, const EntityType& entity, const BaseFunctionSetType& baseFunctionSet ) const
  {
    std::cout << "DiscreteLinearFunctionalInterface::applyLocal()" << std::endl;
    CHECK_INTERFACE_IMPLEMENTATION( asImp().applyLocal( discreteFunctionSpace, entity, baseFunctionSet ) );
    return asImp().applyLocal( discreteFunctionSpace, entity, baseFunctionSet );
  }

protected:

  //! For crtp trick.
  DiscreteLinearFunctionalType& asImp()
  {
      return static_cast< DiscreteLinearFunctionalType& >( *this );
  }

  //! For crtp trick.
  const DiscreteLinearFunctionalType& asImp() const
  {
      return static_cast< const DiscreteLinearFunctionalType& >( *this );
  }

}; // end DiscreteLinearFunctionalInterface

// forward declaration
template< class DiscreteLinearFunctionalDefaultTraitsImp >
class DiscreteLinearFunctionalDefault;

/**
  * \brief      This class is the traits class for the class DiscreteLinearFunctionalDefault.
  *
  * \todo       Doc me, please!
  **/
template< class InducingFunctionImp >
class DiscreteLinearFunctionalDefaultTraits
{
public:

  //! For crtp trick.
  typedef DiscreteLinearFunctionalDefault< DiscreteLinearFunctionalDefaultTraits >
    DiscreteLinearFunctionalType;

  //! Type of the function, which induces the functional.
  typedef InducingFunctionImp
    InducingFunctionType;

}; // end DiscreteLinearFunctionalDefaultTraits


/**
  * \brief      This class is the default implementation of a discrete linear functional.
  *
  *             This class implements the operator() as a gridwalk and calls applyLocal() of the derived class on
  *             each entity.
  *
  * \attention  This class is under construction!
  *
  * \todo       Doc me, please!
  **/
template< class DiscreteLinearFunctionalDefaultTraitsImp >
class DiscreteLinearFunctionalDefault
  : public DiscreteLinearFunctionalInterface< DiscreteLinearFunctionalDefaultTraitsImp >
{
public:

  typedef DiscreteLinearFunctionalDefaultTraitsImp
    Traits;

private:

  typedef DiscreteLinearFunctionalDefault< Traits >
    ThisType;

  typedef DiscreteLinearFunctionalInterface< Traits >
    BaseType;

  //! For crtp trick.
  typedef typename Traits::DiscreteLinearFunctionalType
    DiscreteLinearFunctionalType;

public:

  typedef typename BaseType::LocalDoFVectorType
    LocalDoFVectorType;

  //! Type of the function, which induces the functional.
  typedef typename Traits::InducingFunctionType
    InducingFunctionType;

  //! Type of the inducing functions range field.
  typedef typename InducingFunctionType::RangeFieldType
    RangeFieldType;

  /**
    * \brief    Constructor.
    *
    *           Calls the constructor of the base class.
    **/
  DiscreteLinearFunctionalDefault()
    : BaseType()
  {
    std::cout << "DiscreteLinearFunctionalDefault::DiscreteLinearFunctionalDefault()" << std::endl;
  }

  /**
    * \brief    Destructor (empty).
    **/
  ~DiscreteLinearFunctionalDefault()
  {
    std::cout << "DiscreteLinearFunctionalDefault::~DiscreteLinearFunctionalDefault()" << std::endl;
  }

  /**
    * \brief      This operator represents the application of the functional to a discrete function.
    *
    *             This operator does a grid walk and calls applyLocal() on each entity.
    *
    * \todo       Doc me, please!
    **/
  template< class DiscreteFunctionType >
  const RangeFieldType operator()( const DiscreteFunctionType& discreteFunction ) const
  {
    std::cout << "DiscreteLinearFunctionalDefault::operator()" << std::endl;

    RangeFieldType ret = 0.0;

    // some types we will need
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
      BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType::GridPartType
      GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType
      EntityIteratorType;
    typedef typename EntityIteratorType::Entity
      EntityType;
    typedef typename DiscreteFunctionType::RangeType
      RangeType;
    typedef typename DiscreteFunctionType::LocalFunctionType
      LocalFunctionType;

    // some things we will need
    const DiscreteFunctionSpaceType& discreteFunctionSpace = discreteFunction.space();

    // emulate a gridwalk, do first entity only
    EntityIteratorType entityIterator = discreteFunctionSpace.begin();

    // entity and geometry
    const EntityType& entity = *entityIterator;

    // local function and basefunction set
    const LocalFunctionType& localFunction = discreteFunction.localFunction( entity );
    const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace.baseFunctionSet( entity );

    // local DoF and functional vector
    const LocalDoFVectorType localDoFVector( localFunction );
    const LocalDoFVectorType localFunctionalVector = applyLocal( discreteFunctionSpace, entity, baseFunctionSet );

//      // do gridwalk
//      const EntityIteratorType BehindLastEntityIterator = discreteFunctionSpace.end();
//      for ( EntityIteratorType entityIterator = discreteFunctionSpace.begin();
//            entityIterator != BehindLastEntityIterator;
//            ++entityIterator )
//      {
//        // entity and geometry
//        const EntityType& entity = *entityIterator;

//        // local function and basefunction set
//        const LocalFunctionType& localFunction = discreteFunction.localFunction( entity );
//        const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace.baseFunctionSet( entity );

//        // local DoF and functional vector
//        const LocalDoFVectorType localDoFVector( localFunction );
//        const LocalDoFVectorType localFunctionalVector = applyLocal( discreteFunctionSpace, entity, baseFunctionSet );

//        // compute product
//        ret += localDoFVector * localFunctionalVector;

//      } // done gridwalk

    return ret;
  }

  /**
    * \brief      This function represents the application of the functional to a local basefunctionset.
    *
    *             This function calls the function of the derived class.
    *
    * \todo       Doc me, please!
    **/
  template< class DiscreteFunctionSpaceType, class EntityType, class BaseFunctionSetType >
  const LocalDoFVectorType applyLocal( const DiscreteFunctionSpaceType& discreteFunctionSpace, const EntityType& entity, const BaseFunctionSetType& baseFunctionSet ) const
  {
    std::cout << "DiscreteLinearFunctionalDefault::applyLocal()" << std::endl;
    CHECK_INTERFACE_IMPLEMENTATION( asImp().applyLocal( discreteFunctionSpace, entity, baseFunctionSet ) );
    return asImp().applyLocal( discreteFunctionSpace, entity, baseFunctionSet );
  }

protected:

  //! For crtp trick.
  DiscreteLinearFunctionalType& asImp()
  {
      return static_cast< DiscreteLinearFunctionalType& >( *this );
  }

  //! For crtp trick.
  const DiscreteLinearFunctionalType& asImp() const
  {
      return static_cast< const DiscreteLinearFunctionalType& >( *this );
  }

}; // end class DiscreteLinearFunctionalDefault


// forward declaration
template< class TraitsImp >
class TestLinearFunctional;

/**
  * \brief      This class is the traits class for the class TestLinearFunctional.
  *
  * \attention  This class is only a testclass to toy around with crtp/vritual functions.
  **/
template< class InducingFunctionImp >
class TestLinearFunctionalTraits
{
public:

  //! For crtp trick.
  typedef TestLinearFunctional< TestLinearFunctionalTraits >
    DiscreteLinearFunctionalType;

  //! Type of the function, which induces the functional.
  typedef InducingFunctionImp
    InducingFunctionType;

}; // end DiscreteLinearFunctionalDefaultTraits

/**
  * \brief      This class is only a testclass to toy around with crtp/vritual functions.
  *
  * \attention  This class is only a testclass to toy around with crtp/vritual functions.
  **/
template< class TraitsImp >
class TestLinearFunctional
  : public DiscreteLinearFunctionalDefault< TraitsImp >
{
public:

  typedef TraitsImp
    Traits;

private:

  typedef TestLinearFunctional< TraitsImp >
    ThisType;

  typedef DiscreteLinearFunctionalDefault< TraitsImp >
    BaseType;

public:

  typedef typename Traits::InducingFunctionType
    InducingFunctionType;

  typedef typename BaseType::LocalDoFVectorType
    LocalDoFVectorType;

public:

  //! Constructor
  TestLinearFunctional( const InducingFunctionType& inducingFunction )
    : BaseType(),
      inducingFunction_( inducingFunction )
  {
    std::cout << "TestLinearFunctional::TestLinearFunctional()" << std::endl;
  }

  //! Denstructor
  ~TestLinearFunctional()
  {
    std::cout << "TestLinearFunctional::~TestLinearFunctional()" << std::endl;
  }

  //! Function, that actually does something
  template< class DiscreteFunctionSpaceType, class EntityType, class BaseFunctionSetType >
  const LocalDoFVectorType applyLocal( const DiscreteFunctionSpaceType& discreteFunctionSpace, const EntityType& entity, const BaseFunctionSetType& baseFunctionSet ) const
  {
    std::cout << "TestLinearFunctional::applyLocal()" << std::endl;

    const unsigned size = baseFunctionSet.numBaseFunctions();

    LocalDoFVectorType ret( size );

    return ret;
  }

private:

  const InducingFunctionType& inducingFunction_;

}; // end TestLinearFunctional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_DISCRETELINEARFUNCTIONAL_HH
