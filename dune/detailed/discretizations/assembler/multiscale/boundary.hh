
#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_BOUNDARY_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_BOUNDARY_HH

// system
#include <vector>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/dynmatrix.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Assembler {

namespace Multiscale {

template< class BoundaryGridPartImp,
          class BoundaryInfoImp,
          class AnsatzSpaceImp,
          class TestSpaceImp >
class Boundary
{
public:
  typedef BoundaryGridPartImp BoundaryGridPartType;

  typedef AnsatzSpaceImp AnsatzSpaceType;

  typedef TestSpaceImp TestSpaceType;

  typedef BoundaryInfoImp BoundaryInfoType;

  typedef Boundary< BoundaryGridPartType, BoundaryInfoType, AnsatzSpaceType, TestSpaceType > ThisType;

  Boundary(const BoundaryGridPartType& boundaryGridPart,
           const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo,
           const AnsatzSpaceType& ansatzSpace,
           const TestSpaceType& testSpace)
    : boundaryGridPart_(boundaryGridPart)
    , boundaryInfo_(boundaryInfo)
    , ansatzSpace_(ansatzSpace)
    , testSpace_(testSpace)
  {}

  const BoundaryGridPartType& boundaryGridPart() const
  {
    return boundaryGridPart_;
  }

  const AnsatzSpaceType& ansatzSpace() const
  {
    return ansatzSpace_;
  }

  const TestSpaceType& testSpace() const
  {
    return testSpace_;
  }

  const Dune::shared_ptr< const BoundaryInfoType >boundaryInfo() const
  {
    return boundaryInfo_;
  }

  template< class DirichletAssemblerType, class MatrixBackendType >
  void assemble(const DirichletAssemblerType& dirichletAssembler,
                MatrixBackendType& matrix) const
  {
    // preparations
    typedef typename BoundaryGridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename BoundaryGridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename IntersectionType::EntityPointer EntityPointerType;
    typedef typename AnsatzSpaceType::RangeFieldType RangeFieldType;
    typedef Dune::DynamicMatrix< RangeFieldType > LocalMatrixType;
    // common tmp storage for all entities
    std::vector< unsigned int > numberTmpMatrices = dirichletAssembler.numTmpObjectsRequired();
    std::vector< LocalMatrixType > tmpLocalAssemblerMatrices(numberTmpMatrices[0],
                                                             LocalMatrixType(ansatzSpace_.map().maxLocalSize(),
                                                                             testSpace_.map().maxLocalSize(),
                                                                             RangeFieldType(0.0)));
    std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberTmpMatrices[1],
                                                            LocalMatrixType(ansatzSpace_.map().maxLocalSize(),
                                                                            testSpace_.map().maxLocalSize(),
                                                                            RangeFieldType(0.0)));
    std::vector< std::vector< LocalMatrixType > > tmpLocalMatricesContainer;
    tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
    tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
    // walk the boundary grid part
    for (typename BoundaryGridPartType::template Codim< 0 >::IteratorType entityIt = boundaryGridPart_.template begin< 0 >();
         entityIt != boundaryGridPart_.template end< 0 >();
         ++entityIt) {
      const EntityType& entity = *entityIt;
      // walk the intersections
      for(IntersectionIteratorType interectionIt = boundaryGridPart_.ibegin(entity);
          interectionIt != boundaryGridPart_.iend(entity);
          ++interectionIt ) {
        // get the intersection
        const IntersectionType& intersection = *interectionIt;
        // call the local assembler
        if (boundaryInfo_->dirichlet(intersection)) {
          dirichletAssembler.assembleLocal(intersection,
                                           ansatzSpace_,
                                           testSpace_,
                                           matrix,
                                           tmpLocalMatricesContainer);

        } // call the local assembler
      } // walk the intersections
    } // walk the coupling grid part
  } // void assemble(...) const

private:
  const BoundaryGridPartType& boundaryGridPart_;
  const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo_;
  const AnsatzSpaceType& ansatzSpace_;
  const TestSpaceType& testSpace_;
}; // class Boundary

} // namespace Multiscale

} // namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_BOUNDARY_HH
