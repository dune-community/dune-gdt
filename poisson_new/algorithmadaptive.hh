#ifndef POISSON_ADAPTIVE_HH
#define POISSON_ADAPTIVE_HH


#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/lagrangespace/adaptmanager.hh>

#include "algorithm.hh"
#include "estimator.hh"

// Algorithm
// --------- 
template <class GridImp, int polynomialOrder> 
class AdaptiveAlgorithm :
  public Algorithm<GridImp, 
                   polynomialOrder,
                   Dune::AdaptiveLeafGridPart< GridImp >  
                  > /*@LST0S@*/
{
  typedef Algorithm<GridImp, polynomialOrder, Dune::AdaptiveLeafGridPart< GridImp > > BaseType;
public:
  //------ HGridType ----------------------------------------------------------
  typedef typename BaseType :: HGridType HGridType;

  //------ GridPartType ----------------------------------------------------------
  typedef typename BaseType :: GridPartType GridPartType;

  //! to be revised 
  typedef typename BaseType :: ProblemType ProblemType;

  //---- DiscreteFunction ----------------------------------------------------
  typedef typename BaseType :: DiscreteFunctionType  DiscreteFunctionType;

  //---- Local Restriction and Prolongation Operator -------------------------
  typedef Dune::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
  //---- Adaptation Manager --------------------------------------------------
  typedef Dune::AdaptationManager< HGridType, RestrictionProlongationType > AdaptationManagerType;

  using BaseType :: problem_;
public:
  //! constructor
  AdaptiveAlgorithm(HGridType & grid, const ProblemType& problem) :          /*@LST0E@*/
    BaseType( grid, problem )
  {
    if( ! GridPartType :: conforming ) 
    {
      std::cerr << "Non-conforming refinement not supported in this tutorial!" <<std::endl;
      std::cerr << "Choose a grid with conforming refinement, such as AlbertaGrid. " << std::endl;
      abort();
    }
  }

  //! estimate and mark solution 
  bool estimateAndMark(const DiscreteFunctionType& solution, const double tolerance) const 
  {
    std::cout << "Got tolerance = " << tolerance << std::endl;
    Estimator< DiscreteFunctionType > estimator( solution );
    // estimate and mark 
    return estimator.estimateAndMark( problem_, tolerance ); 
  }
};

#endif // POISSONADAPTIVE_HH
