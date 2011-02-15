#ifndef DUNE_FEM_FUNCTIONALS_FUNCTIONTOOLS_HH
#define DUNE_FEM_FUNCTIONALS_FUNCTIONTOOLS_HH

// dune grid includes
#include <dune/fem/io/file/vtkio.hh>

namespace Dune
{

namespace FemTools{

  /**
    * \brief  Function to set a discrete function to a constant value.
    *
    *         There should be an easier way to do this in Fem, but i can't come up with one at the moment.
    **/
  template< class DiscreteFunctionType >
  void setDiscreteFunctionToScalarValue( DiscreteFunctionType& discreteFunction, const double& scalar )
  {
    discreteFunction.clear();

    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
    DofIteratorType it = discreteFunction.dbegin();
    DofIteratorType itEnd = discreteFunction.dend();
    for ( ; it != itEnd; ++it )
        *it = scalar;
    return;
  }

  /**
    * \brief  Function to write a discrete function to .vtk for visualizytion purpose.
    **/
  template< class DiscreteFunctionType >
  void writeDiscreteFunctionToVTK( const DiscreteFunctionType& discreteFunction, const std::string filename = "discrete_function_written_by_fem_tools" )
  {
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType::GridPartType
      GridPartType;
    const GridPartType& gridPart = discreteFunction.space().gridPart();

    typedef Dune::VTKIO< GridPartType >
      VTKWriterType;
    VTKWriterType vtkWriter( gridPart );

//    vtkWriter.addVectorVertexData( computedSolutions.discreteVelocity() );
    vtkWriter.addVertexData( discreteFunction );
    vtkWriter.write( filename );
    vtkWriter.clear();

  }

} // end namespace FemTools

} // end namespace Dune




#endif // DUNE_FEM_FUNCTIONALS_FUNCTIONTOOLS_HH
