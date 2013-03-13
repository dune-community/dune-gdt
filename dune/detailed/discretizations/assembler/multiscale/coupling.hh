#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_COUPLING_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_COUPLING_HH

#include <vector>
#include <memory>

#include <dune/common/dynmatrix.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Assembler {

namespace Multiscale {

namespace Coupling {

template< class CouplingGridPartImp,
          class InnerAnsatzSpaceImp,
          class InnerTestSpaceImp,
          class OuterAnsatzSpaceImp,
          class OuterTestSpaceImp >
class Primal
{
public:
  typedef CouplingGridPartImp CouplingGridPartType;

  typedef InnerAnsatzSpaceImp InnerAnsatzSpaceType;

  typedef InnerTestSpaceImp InnerTestSpaceType;

  typedef OuterAnsatzSpaceImp OuterAnsatzSpaceType;

  typedef OuterTestSpaceImp OuterTestSpaceType;

  typedef Primal< CouplingGridPartType, InnerAnsatzSpaceType, InnerTestSpaceType, OuterAnsatzSpaceType, OuterTestSpaceType > ThisType;

private:
  typedef typename CouplingGridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename InnerAnsatzSpaceType::RangeFieldType RangeFieldType;
  typedef Dune::DynamicMatrix< RangeFieldType > LocalMatrixType;
  typedef std::vector< std::vector< LocalMatrixType > > LocalMatricesContainerType;

  class LocalMatrixAssemblerApplication
  {
  public:
    virtual ~LocalMatrixAssemblerApplication(){}

    virtual void apply(const IntersectionType& /*_intersection*/,
                       const InnerAnsatzSpaceType& /*_innerAnsatzSpace*/,
                       const InnerTestSpaceType& /*_innerTestSpace*/,
                       const OuterAnsatzSpaceType& /*_outerAnsatzSpace*/,
                       const OuterTestSpaceType& /*_outerTestSpace*/,
                       LocalMatricesContainerType& /*_localMatricesContainer*/) const = 0;

    virtual std::vector< unsigned int > numTmpObjectsRequired() const = 0;
  }; // class LocalOperatorApplication

  template< class LocalMatrixAssemblerType, class MatrixType >
  class LocalMatrixAssemblerApplicationWrapper
    : public LocalMatrixAssemblerApplication
  {
  public:
    LocalMatrixAssemblerApplicationWrapper(const std::shared_ptr< const LocalMatrixAssemblerType > _localMatrixAssembler,
                                           std::shared_ptr< MatrixType > _innerInnerMatrix,
                                           std::shared_ptr< MatrixType > _outerOuterMatrix,
                                           std::shared_ptr< MatrixType > _innerOuterMatrix,
                                           std::shared_ptr< MatrixType > _outerInnerMatrix)
      : localMatrixAssembler_(_localMatrixAssembler)
      , innerInnerMatrix_(_innerInnerMatrix)
      , outerOuterMatrix_(_outerOuterMatrix)
      , innerOuterMatrix_(_innerOuterMatrix)
      , outerInnerMatrix_(_outerInnerMatrix)
    {}

    virtual void apply(const IntersectionType& _intersection,
                       const InnerAnsatzSpaceType& _innerAnsatzSpace,
                       const InnerTestSpaceType& _innerTestSpace,
                       const OuterAnsatzSpaceType& _outerAnsatzSpace,
                       const OuterTestSpaceType& _outerTestSpace,
                       LocalMatricesContainerType& _localMatricesContainer) const
    {
      localMatrixAssembler_->assembleLocal(_intersection,
                                           _innerAnsatzSpace,
                                           _innerTestSpace,
                                           _outerAnsatzSpace,
                                           _outerTestSpace,
                                           *innerInnerMatrix_,
                                           *outerOuterMatrix_,
                                           *innerOuterMatrix_,
                                           *outerInnerMatrix_,
                                           _localMatricesContainer);
    } // virtual void applyLocal(...) const

    virtual std::vector< unsigned int > numTmpObjectsRequired() const
    {
      return localMatrixAssembler_->numTmpObjectsRequired();
    } // virtual std::vector< unsigned int > numTmpObjectsRequired() const

  private:
    const std::shared_ptr< const LocalMatrixAssemblerType > localMatrixAssembler_;
    std::shared_ptr< MatrixType > innerInnerMatrix_;
    std::shared_ptr< MatrixType > outerOuterMatrix_;
    std::shared_ptr< MatrixType > innerOuterMatrix_;
    std::shared_ptr< MatrixType > outerInnerMatrix_;
  }; // class LocalMatrixAssemblerApplicationWrapper

public:
  Primal(const CouplingGridPartType& couplingGridPart,
         const InnerAnsatzSpaceType& innerAnsatzSpace,
         const InnerTestSpaceType& innerTestSpace,
         const OuterAnsatzSpaceType& outerAnsatzSpace,
         const OuterTestSpaceType& outerTestSpace)
    : couplingGridPart_(couplingGridPart)
    , innerAnsatzSpace_(innerAnsatzSpace)
    , innerTestSpace_(innerTestSpace)
    , outerAnsatzSpace_(outerAnsatzSpace)
    , outerTestSpace_(outerTestSpace)
  {}

  ~Primal()
  {
    for (auto& localMatrixAssembler: localMatrixAssemblers_)
      delete localMatrixAssembler;
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

  template< class LocalMatrixAssemblerType, class MatrixType >
  void addLocalMatrixAssembler(const std::shared_ptr< const LocalMatrixAssemblerType > _localMatrixAssembler,
                               std::shared_ptr< MatrixType > _innerInnerMatrix,
                               std::shared_ptr< MatrixType > _outerOuterMatrix,
                               std::shared_ptr< MatrixType > _innerOuterMatrix,
                               std::shared_ptr< MatrixType > _outerInnerMatrix)
  {
    typedef LocalMatrixAssemblerApplicationWrapper< LocalMatrixAssemblerType, MatrixType > WrapperType;
    WrapperType* wrapper = new WrapperType(_localMatrixAssembler,
                                           _innerInnerMatrix,
                                           _outerOuterMatrix,
                                           _innerOuterMatrix,
                                           _outerInnerMatrix);
    localMatrixAssemblers_.push_back(wrapper);
  }

//  template< class LocalAssemblerType, class MatrixBackendType >
  void assemble/*Matrices*/(/*const LocalAssemblerType& localAssembler,
                        MatrixBackendType& innerInnerMatrix,
                        MatrixBackendType& innerOuterMatrix,
                        MatrixBackendType& outerInnerMatrix,
                        MatrixBackendType& outerOuterMatrix*/) const
  {
    // preparations
    typedef typename CouplingGridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename IntersectionType::EntityPointer EntityPointerType;
    // common tmp storage for all entities
    // * for the matrix assemblers
    std::vector< unsigned int > numberOfTmpMatricesNeeded(2, 0);
    for (unsigned int ii = 0; ii < localMatrixAssemblers_.size(); ++ii) {
      const std::vector< unsigned int > tmp = localMatrixAssemblers_[ii]->numTmpObjectsRequired();
      numberOfTmpMatricesNeeded[0] = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
      numberOfTmpMatricesNeeded[1] = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
    }
    std::vector< LocalMatrixType > tmpLocalAssemblerMatrices(numberOfTmpMatricesNeeded[0],
                                                             LocalMatrixType(std::max(innerAnsatzSpace_.map().maxLocalSize(), outerAnsatzSpace_.map().maxLocalSize()),
                                                                             std::max(innerTestSpace_.map().maxLocalSize(), outerTestSpace_.map().maxLocalSize()),
                                                                             RangeFieldType(0.0)));
    std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberOfTmpMatricesNeeded[1],
                                                            LocalMatrixType(std::max(innerAnsatzSpace_.map().maxLocalSize(), outerAnsatzSpace_.map().maxLocalSize()),
                                                                            std::max(innerTestSpace_.map().maxLocalSize(), outerTestSpace_.map().maxLocalSize()),
                                                                            RangeFieldType(0.0)));
    std::vector< std::vector< LocalMatrixType > > tmpLocalMatricesContainer;
    tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
    tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
//    // common tmp storage for all entities
//    std::vector< unsigned int > numberTmpMatrices = localAssembler.numTmpObjectsRequired();
//    std::vector< LocalMatrixType > tmpLocalAssemblerMatrices(numberTmpMatrices[0],
//                                                             LocalMatrixType(std::max(innerAnsatzSpace_.map().maxLocalSize(), outerAnsatzSpace_.map().maxLocalSize()),
//                                                                             std::max(innerTestSpace_.map().maxLocalSize(), outerTestSpace_.map().maxLocalSize()),
//                                                                             RangeFieldType(0.0)));
//    std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberTmpMatrices[1],
//                                                            LocalMatrixType(std::max(innerAnsatzSpace_.map().maxLocalSize(), outerAnsatzSpace_.map().maxLocalSize()),
//                                                                            std::max(innerTestSpace_.map().maxLocalSize(), outerTestSpace_.map().maxLocalSize()),
//                                                                            RangeFieldType(0.0)));
//    std::vector< std::vector< LocalMatrixType > > tmpLocalMatricesContainer;
//    tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
//    tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
    // walk the coupling grid part
    for (typename CouplingGridPartType::template Codim< 0 >::IteratorType entityIt = couplingGridPart_.template begin< 0 >();
         entityIt != couplingGridPart_.template end< 0 >();
         ++entityIt) {
      const EntityType& insideEntity = *entityIt;
      // walk the intersections
      for(IntersectionIteratorType interectionIt = couplingGridPart_.ibegin(insideEntity);
          interectionIt != couplingGridPart_.iend(insideEntity);
          ++interectionIt ) {
        // with a coupling grid part we can be sure to only get inner intersection (so we only need the assert, not an if)
        const IntersectionType& intersection = *interectionIt;
        assert(intersection.neighbor() && !intersection.boundary());
        // call the local assembler
        for (unsigned int ii = 0; ii < localMatrixAssemblers_.size(); ++ii)
          localMatrixAssemblers_[ii]->apply(intersection,
                                            innerAnsatzSpace_,
                                            innerTestSpace_,
                                            outerAnsatzSpace_,
                                            outerTestSpace_,
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
  std::vector< LocalMatrixAssemblerApplication* > localMatrixAssemblers_;
}; // class Primal

} // namespace Coupling

} // namespace Multiscale

} // namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_COUPLING_HH
