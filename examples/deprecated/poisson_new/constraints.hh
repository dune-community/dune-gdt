#ifndef DUNE_CONSTRAINTS
#define DUNE_CONSTRAINTS

#include "functional.hh"

/* 
 * @class A class representing a single constraint
 *
 * A single constraint knows the constraint defined by a functional, the row number
 * an assembler has to delete in the system matrix. The constraint can be either
 * applied to a discrete function or to a local basis function.
 * 
 * @todo The template argument should be a DiscreteFunctionalType by now. Change it, if needed.
 */
template< class FunctionalType >
class SingleConstraint
{
  public:
  typedef typename FunctionalType::Codim
    Codim;
  typedef typename FunctionalType::DiscreteFunctionSpace
    DiscreteFunctionSpace;

  private:
  typedef typename DiscreteFunctionSpace::BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename DiscreteFunctionSpace::DiscreteFunction
    DiscreteFunctionType;
  typedef typename DiscreteFunctionSpace::Traits
    DFSTraits;

  typedef typename DFSTraits::GridPartType
    GridPartType;


  typedef typename GridPartType::GridType
    GridType;
  typedef typename GridPartType::IndexSetType
    IndexSetType;
  typedef typename GridPartType::template Codim< 0 >::IteratorType
    LeafIterator;
  typedef typename GridPartType::IntersectionIterator
    IntersectionIterator;
  typedef typename GridPartType::IndexSet
    LeafIndexSet;
  typedef typename LeafIterator::Entity
    EntityType;

  public:


  /*
   * @brief constructor
   *
   * @ f The Functional describing the condition
   * @ row The row which has to be deleted by the assembler. The default value means that no row is deleted
   */
  SingleConstraint( FunctionalType& f, int row=-1 ):f_( f ), row_( row )
  {

  }

  /*
   * @brief Returns the functional for this constraint
   *
   */
  FunctionalType getFunctional()
  {
    return f_;
  }


  /*
   * @brief Returns the row which has to be deleted by the assembler
   *
   * The value -1 means that no row is deleted
   *
   */
  int getRow()
  {
    return row_;
  }

  /*
   * @brief Applys the discrete function to the constraint defined by the functional returned by getFunctional() 
   *
   * We want to solve the following problem: Find a discrete function u_h which fulfills the condition 
   * f(u_h)=0 for the functional f given by getFunctional().
   *
   * @param u_h the discrete function
   * 
   * @todo implement apply(FunctionType& u), if needed 
   */
  void apply( DiscreteFunctionType& u_h )
  {

  }

  /*
   * @brief Applys the local base function bf on entity en to the constraint defined 
   * by the functional returned by getFunctional() 
   *
   * @param en the entity
   * @param bf the local base function
   *
   * We want to solve the following problem: Find a local base function bf which fulfills the condition 
   * f(bf)=0 on entity en for the functional f given by getFunctional().
   *
   * @todo How many variables of BaseFunctionSetType are needed?
   * @todo How to deal with the sparse vector? What requirements are needed?
   *
   */
  void applyLocal( EntityType& en, BaseFunctionSetType& bf )
  {

  }

  private:

  int row_;
  FunctionalType f_;

};


/*
 * @class A class all constraints classes have to be derived from
 *
 */
template< class FunctionalType >
class ConstraintsInterface
{
  public:
  typedef typename FunctionalType::Codim
    Codim;
  typedef typename FunctionalType::DiscreteFunctionSpace
    DiscreteFunctionSpace;
  typedef  SingleConstraint< FunctionalType >
    SingleConstraintType;

  private:
  typedef typename DiscreteFunctionSpace::BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename DiscreteFunctionSpace::DiscreteFunction
    DiscreteFunctionType;
  typedef typename DiscreteFunctionSpace::Traits
    DFSTraits;

  typedef typename DFSTraits::GridPartType
    GridPartType;


  typedef typename GridPartType::GridType
    GridType;
  typedef typename GridPartType::IndexSetType
    IndexSetType;
  typedef typename GridPartType::template Codim< 0 >::IteratorType
    LeafIterator;
  typedef typename GridPartType::IntersectionIterator
    IntersectionIterator;
  typedef typename GridPartType::IndexSet
    LeafIndexSet;

  typedef typename LeafIterator::Entity
    EntityType;


  public:

    /*
     * @brief Constructor
     */
    ConstraintsInterface()
    {
    }

    /*
     * @brief returns number of constraints
     */
    virtual std::size_t size() = 0;



    /*
     * @brief Applys the discrete function to the constraint defined by the functional returned by getFunctional() 
     *
     * We want to solve the following problem: Find a discrete function u_h which fulfills the condition 
     * f(u_h)=0 for the functional f given by getFunctional().
     * 
     * Normally this is only a loop over all constraints applied to the discrete function.
     *
     * @param u_h the discrete function
     */
    virtual void apply( DiscreteFunctionType u_h ) = 0;
  
    /*
     * @brief apply all constraints on Entity en and fills sparse matrix 
     * applys the local base function bf on entity en to the constraint defined 
     * by the functional returned by getFunctional() 
     *
     * @param en the entity
     * @param bf the local base function
     *
     * We want to solve the following problem: Find a local base function bf which fulfills the condition 
     * f(bf)=0 on entity en for the functional f given by getFunctional().
     * 
     * Normally this is only a loop over all constraints applied to the basis function.
     *
     * @todo How many variables of BaseFunctionSetType are needed?
     * @todo How to deal with the sparse matrix? What requirements are needed?
     */
    virtual void applyLocal( EntityType& en, BaseFunctionSetType bf ) = 0;

    /*
     * @brief Returns a single constraint 
     *
     * @param i (Positive) number of constraint
     *
     * This operator returns a single constraint. 
     */
    virtual SingleConstraintType operator[]( unsigned int i ) = 0;

};



/*
 * @class Default implementation of the class ConstraintsInterface
 *
 * @todo add methods
 */
template< class FunctionalType >
class ConstraintsDefault
  : public ConstraintsInterface< FunctionalType >
{
public:
  //SingleConstraintType operator[](unsigned int i){
  //  return nil;
  //}

private:

};


/*
 * @class Constraints class implementing Dirichlet constraints on boundary 
 * 
 * @todo add methods
 */
template< class FunctionalType >
class DirichletBoundaryConstraints
  : public ConstraintsDefault< FunctionalType >
{ 
public:


  std::size_t size()
  {
    //do a grid walk and count boundary faces
    return 1;
  }

  //In general we are not really keen on using this method here:
  //TODO Do your own implementation here without using the [] operator
  void apply( DiscreteFunctionType& u_h )
  {
    //TODO remove it
    for ( int i = 0; i < size(); i++ )
      ( *this )[i].getFunctional().apply( u_h );
  }


  //In general we are not really keen on using this method here:
  //TODO Do your own implementation here without using the [] operator
  void applyLocal( EntityType& en, BaseFunctionSetType& bf )
  {
    //TODO remove it
    for ( int i = 0; i < size(); i++ )
      ( *this )[i].getFunctional().applyLocal( bf );
  }



  //In general we are not really keen on using this method here
  SingleConstraintType operator[]( unsigned int i )
  {
    //Do a grid walk here and take the i-th boundary faces
    //return a single constraint
    
    //create the right Functional deleting the correct dof
    FunctionalType functional;
    int dummyline = 1;
    SingleConstraint constraint( functional, dummyline );
    return constraint;
  }

private:

};


/*
 * @class Constraints class implementing Neumann constraints on boundary
 *
 * @todo add methods
 */
template< class FunctionalType >
class NeumannBoundaryConstraints
  : public ConstraintsDefault< FunctionalType >{

};



/*
 * @class Constraint class where constraints are stored in a vector
 *
 * @todo implement methods
 */
template< class FunctionalType >
class VectorConstraints
  : public ConstraintsDefault< FunctionalType >
{
  public:
  typedef typename FunctionalType::Codim
    Codim;
  typedef typename FunctionalType::DiscreteFunctionSpace
    DiscreteFunctionSpace;
  typedef SingleConstraint< FunctionalType >
    SingleConstraintType;

  private:
  typedef typename DiscreteFunctionSpace::BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename DiscreteFunctionSpace::DiscreteFunction
    DiscreteFunctionType;
  typedef typename DiscreteFunctionSpace::Traits
    DFSTraits;

  typedef typename DFSTraits::GridPartType
    GridPartType;

  typedef typename GridPartType::GridType
    GridType;
  typedef typename GridPartType::IndexSetType
    IndexSetType;
  typedef typename GridPartType::template Codim< 0 >::IteratorType
    LeafIterator;
  typedef typename GridPartType::IntersectionIterator
    IntersectionIterator;
  typedef typename GridPartType::IndexSet
    LeafIndexSet;

  typedef typename LeafIterator::Entity
    EntityType;


  public:

    /*
     * @brief Constructor
     *
     * @param constraint a vector storing the constraints
     */
    VectorConstraints(std::vector<SingleConstraintType> constraints)
    {
    } 

    /*
     * @brief returns number of constraints
     */
    std::size_t size()
    {
      return constraints_.size();
    }

    /*
     * @brief Applys the discrete function to the constraint defined by the functional returned by getFunctional() 
     *
     * We want to solve the following problem: Find a discrete function u_h which fulfills the condition 
     * f(u_h)=0 for the functional f given by getFunctional().
     * 
     * Normally this is only a loop over all constraints applied to the discrete function.
     *
     * @param u_h the discrete function
     */
    void apply( DiscreteFunctionType& u_h )
    {
      for (int i = 0; i < size(); i++ )
        constraints_[i].getFunctional().apply( u_h );
    }

    /*
     * @brief apply all constraints on Entity en and fills sparse matrix 
     * applys the local base function bf on entity en to the constraint defined 
     * by the functional returned by getFunctional() 
     *
     * @param en the entity
     * @param bf the local base function
     *
     * We want to solve the following problem: Find a local base function bf which fulfills the condition 
     * f(bf)=0 on entity en for the functional f given by getFunctional().
     * 
     * Normally this is only a loop over all constraints applied to the basis function.
     *
     * @todo How many variables of BaseFunctionSetType are needed?
     * @todo How to deal with the sparse matrix? What requirements are needed?
     */
    void applyLocal( EntityType& en, BaseFunctionSetType& bf )
    {
      for ( int i = 0; i < size(); i++ )
        constraints_[i].getFunctional().applyLocal( bf );
    }


    /*
    SparsityPattern localSparsityPattern(EntityType& en){

      typedef std::vector<std::set<int> >                    SparsityPattern;

      Dune::GeometryType gt = en.type();

      const Dune::template GenericReferenceElement<ctype,dim> &ref =
         Dune::GenericReferenceElements<ctype,dim>::general(gt);

      // traverse all codim-1-entities of the current element and store all
      // pairs of vertices in sparsityPattern_
      const IntersectionIterator isend = en.iend(en);
      for (IntersectionIterator is = en.ibegin(en) ; is != isend ; ++is)
      {
        int vertexsize = ref.size(is->indexInInside(),1,dim);
        for (int i=0; i < vertexsize; i++)
        {
          int indexi = set.subIndex(en,ref.subEntity(is->indexInInside(),1,i,dim),dim);
          for (int j=0; j < vertexsize; j++)
          {
            int indexj = set.subIndex(en,ref.subEntity(is->indexInInside(),1,j,dim),dim);
            sparsityPattern_[indexi].insert(indexj);
          }
        }
       }

    }*/

    /*
     * @brief Returns a single constraint 
     *
     * @param i (Positive) number of constraint
     *
     * This operator returns a single constraint. 
     */
    SingleConstraintType operator[]( unsigned int i )
    {
      return constraints_[i];
    }

  private:

    /*
     * @brief checks whether all nonnegative numbers are unique
     */
    bool checkConstraintMapUnique()
    {
      return true;
    }

    std::vector< SingleConstraintType > constraints_;


};

#endif
