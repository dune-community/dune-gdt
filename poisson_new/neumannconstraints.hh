#ifndef DUNE_NEUMANNCONSTRAINTS_HH
#define DUNE_NEUMANNCONSTRAINTS_HH


namespace Dune {


template <class DiscreteFunctionSpace, class ProblemType>
class NeumannConstraints
{
public:
  typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

  //! type of grid partition
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  //! type of grid
  typedef typename DiscreteFunctionSpaceType::GridType GridType;

  NeumannConstraints(const DiscreteFunctionSpaceType& space, const ProblemType& problem)
    : space_(space)
    , problem_(problem)
  {
  }

  template <class LinearOperator, class DiscreteFunctionType>
  void apply(LinearOperator& linearOperator, DiscreteFunctionType& rhs, DiscreteFunctionType& solution) const
  {
    typedef typename ProblemType::ExactSolutionType ExactSolutionType;
    //---- Adapter for exact solution ------------------------------------------
    typedef Dune::DiscreteFunctionAdapter<ExactSolutionType, GridPartType> GridExactSolutionType;

    // create adapter (for visualization with grape)
    GridExactSolutionType ugrid(
        "exact solution", problem_.exactSolution(), space_.gridPart(), DiscreteFunctionSpaceType::polynomialOrder + 1);

    // if no dirichlet boundary is set then assume Neumann boundary
    clearMeanValue(rhs);
    setExactMeanValue(ugrid, solution);
  }

  //! substract the mean value from the right hand side
  template <class DiscreteFunctionType>
  void clearMeanValue(DiscreteFunctionType& rhs) const
  {
    /*@LST0E@*/
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;

    typename DiscreteFunctionSpaceType::RangeFieldType meanValue(0);

    const DofIteratorType dend = rhs.dend();
    for (DofIteratorType dit = rhs.dbegin(); dit != dend; ++dit)
      meanValue += *dit;
    meanValue /= rhs.size();

    std::cout << "Substracting mean value from right hand side: " << meanValue << std::endl;
    for (DofIteratorType dit = rhs.dbegin(); dit != dend; ++dit)
      *dit -= meanValue;
    /*@LST0S@*/
  }

  //! set the exact mean value to the solution
  template <class GridFunctionType, class DiscreteFunctionType>
  void setExactMeanValue(const GridFunctionType& ugrid, DiscreteFunctionType& solution) const
  {
    /*@LST0E@*/
    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;

    typename FunctionSpaceType::RangeType meanValue(0);
    typename FunctionSpaceType::RangeFieldType volume(0);

    const IteratorType end = space_.end();
    for (IteratorType it = space_.begin(); it != end; ++it) {
      const typename IteratorType::Entity& entity               = *it;
      const typename IteratorType::Entity::Geometry& geometry   = entity.geometry();
      const typename GridFunctionType::LocalFunctionType ulocal = ugrid.localFunction(entity);

      //! type of quadrature to be used
      typedef CachingQuadrature<GridPartType, 0> QuadratureType;
      QuadratureType quadrature(entity, space_.order() + 1);
      const int numQuadraturePoints = quadrature.nop();
      for (int qp = 0; qp < numQuadraturePoints; ++qp) {
        const typename FunctionSpaceType::RangeFieldType factor =
            quadrature.weight(qp) * geometry.integrationElement(quadrature.point(qp));
        volume += factor;

        typename FunctionSpaceType::RangeType value;
        ulocal.evaluate(quadrature[qp], value);
        meanValue.axpy(factor, value);
      }
    }
    meanValue /= volume;

    std::cout << "Setting exact mean value: " << meanValue[0] << std::endl;
    const DofIteratorType dend = solution.dend();
    for (DofIteratorType dit = solution.dbegin(); dit != dend; ++dit)
      *dit = meanValue[0];
    /*@LST0S@*/
  }

protected:
  //! pointer to slave dofs
  const DiscreteFunctionSpaceType& space_;
  const ProblemType& problem_;
};

} // end namespace
#endif
