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
    * \brief      This class is the interfac for all discrete linear functionals.
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
    RangeFieldType operator()( const DiscreteFunctionType& discreteFunction ) const
    {
      std::cout << "DiscreteLinearFunctionalInterface::operator()" << std::endl;
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator()( discreteFunction ) );
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

  public:

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
      *             This operator does a grid walk and calls applyLocal() of the derived class on each entity.
      *
      * \todo       Doc me, please!
      **/
    template< class DiscreteFunctionType >
    RangeFieldType operator()( const DiscreteFunctionType& discreteFunction ) const
    {
      RangeFieldType ret = 0.0;

      std::cout << "DiscreteLinearFunctionalDefault::operator()" << std::endl;

      return ret;
    }

  }; // end class DiscreteLinearFunctionalDefault

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_DISCRETELINEARFUNCTIONAL_HH
