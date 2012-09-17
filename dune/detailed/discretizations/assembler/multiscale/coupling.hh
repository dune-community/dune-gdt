
#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_COUPLING_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_COUPLING_HH

// system
#include <vector>

// dune-common
#include <dune/common/dynmatrix.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Assembler {

namespace Multiscale {

namespace Coupling {

template <class CouplingGridPartImp, class InnerAnsatzSpaceImp, class InnerTestSpaceImp, class OuterAnsatzSpaceImp,
          class OuterTestSpaceImp>
class Primal
{
public:
  typedef CouplingGridPartImp CouplingGridPartType;

  typedef InnerAnsatzSpaceImp InnerAnsatzSpaceType;

  typedef InnerTestSpaceImp InnerTestSpaceType;

  typedef OuterAnsatzSpaceImp OuterAnsatzSpaceType;

  typedef OuterTestSpaceImp OuterTestSpaceType;

  typedef Primal<CouplingGridPartType, InnerAnsatzSpaceType, InnerTestSpaceType, OuterAnsatzSpaceType,
                 OuterTestSpaceType> ThisType;

  Primal(const CouplingGridPartType& couplingGridPart, const InnerAnsatzSpaceType& innerAnsatzSpace,
         const InnerTestSpaceType& innerTestSpace, const OuterAnsatzSpaceType& outerAnsatzSpace,
         const OuterTestSpaceType& outerTestSpace)
    : couplingGridPart_(couplingGridPart)
    , innerAnsatzSpace_(innerAnsatzSpace)
    , innerTestSpace_(innerTestSpace)
    , outerAnsatzSpace_(outerAnsatzSpace)
    , outerTestSpace_(outerTestSpace)
  {
  }

  const CouplingGridPartType& couplingGridPart() const
  {
    return couplingGridPart_;
  }

  const InnerAnsatzSpaceType& innerAnsatzSpace() const
  {
    return innerAnsatzSpace_;
  }

  const InnerTestSpaceType& innerTestSpace() const
  {
    return innerTestSpace_;
  }

  const OuterAnsatzSpaceType& outerAnsatzSpace() const
  {
    return outerAnsatzSpace_;
  }

  const OuterTestSpaceType& outerTestSpace() const
  {
    return outerTestSpace_;
  }

  template <class LocalAssemblerType, class MatrixBackendType>
  void assembleMatrices(const LocalAssemblerType& localAssembler, MatrixBackendType& innerInnerMatrix,
                        MatrixBackendType& innerOuterMatrix, MatrixBackendType& outerInnerMatrix,
                        MatrixBackendType& outerOuterMatrix) const
  {
    // preparations
    typedef typename CouplingGridPartType::template Codim<0>::EntityType EntityType;
    typedef typename CouplingGridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename IntersectionType::EntityPointer EntityPointerType;
    typedef typename InnerAnsatzSpaceType::RangeFieldType RangeFieldType;
    typedef Dune::DynamicMatrix<RangeFieldType> LocalMatrixType;
    // common tmp storage for all entities
    std::vector<unsigned int> numberTmpMatrices = localAssembler.numTmpObjectsRequired();
    std::vector<LocalMatrixType> tmpLocalAssemblerMatrices(
        numberTmpMatrices[0],
        LocalMatrixType(std::max(innerAnsatzSpace_.map().maxLocalSize(), outerAnsatzSpace_.map().maxLocalSize()),
                        std::max(innerTestSpace_.map().maxLocalSize(), outerTestSpace_.map().maxLocalSize()),
                        RangeFieldType(0.0)));
    std::vector<LocalMatrixType> tmpLocalOperatorMatrices(
        numberTmpMatrices[1],
        LocalMatrixType(std::max(innerAnsatzSpace_.map().maxLocalSize(), outerAnsatzSpace_.map().maxLocalSize()),
                        std::max(innerTestSpace_.map().maxLocalSize(), outerTestSpace_.map().maxLocalSize()),
                        RangeFieldType(0.0)));
    std::vector<std::vector<LocalMatrixType>> tmpLocalMatricesContainer;
    tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
    tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
    // walk the coupling grid part
    for (typename CouplingGridPartType::template Codim<0>::IteratorType entityIt =
             couplingGridPart_.template begin<0>();
         entityIt != couplingGridPart_.template end<0>();
         ++entityIt) {
      const EntityType& insideEntity = *entityIt;
      // walk the intersections
      for (IntersectionIteratorType interectionIt = couplingGridPart_.ibegin(insideEntity);
           interectionIt != couplingGridPart_.iend(insideEntity);
           ++interectionIt) {
        // with a coupling grid part we can be sure to only get inner intersection (so we only need the assert, not an
        // if)
        const IntersectionType& intersection = *interectionIt;
        assert(intersection.neighbor() && !intersection.boundary());
        // call the local assembler
        localAssembler.assembleLocal(intersection,
                                     innerAnsatzSpace_,
                                     innerTestSpace_,
                                     outerAnsatzSpace_,
                                     outerTestSpace_,
                                     innerInnerMatrix,
                                     outerOuterMatrix,
                                     innerOuterMatrix,
                                     outerInnerMatrix,
                                     tmpLocalMatricesContainer);
      } // walk the intersections
    } // walk the coupling grid part
  } // void assembleMatrices() const

private:
  const CouplingGridPartType& couplingGridPart_;
  const InnerAnsatzSpaceType& innerAnsatzSpace_;
  const InnerTestSpaceType& innerTestSpace_;
  const OuterAnsatzSpaceType& outerAnsatzSpace_;
  const OuterTestSpaceType& outerTestSpace_;
}; // class Primal

} // namespace Coupling

} // namespace Multiscale

} // namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_COUPLING_HH
