
#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_BOUNDARY_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_BOUNDARY_HH

#include <vector>
#include <memory>

#include <dune/common/dynmatrix.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Assembler {

namespace Multiscale {

template <class BoundaryGridPartImp, class BoundaryInfoImp, class AnsatzSpaceImp, class TestSpaceImp>
class Boundary
{
public:
  typedef BoundaryGridPartImp BoundaryGridPartType;

  typedef AnsatzSpaceImp AnsatzSpaceType;

  typedef TestSpaceImp TestSpaceType;

  typedef BoundaryInfoImp BoundaryInfoType;

  typedef Boundary<BoundaryGridPartType, BoundaryInfoType, AnsatzSpaceType, TestSpaceType> ThisType;

private:
  typedef typename BoundaryGridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename AnsatzSpaceType::RangeFieldType RangeFieldType;
  typedef Dune::DynamicMatrix<RangeFieldType> LocalMatrixType;
  typedef std::vector<std::vector<LocalMatrixType>> LocalMatricesContainerType;

  class LocalMatrixAssemblerApplication
  {
  public:
    virtual ~LocalMatrixAssemblerApplication()
    {
    }

    virtual void apply(const AnsatzSpaceType& /*ansatzSpace_*/, const TestSpaceType& /*testSpace_*/,
                       const IntersectionType& /*intersection*/,
                       LocalMatricesContainerType& /*_localMatricesContainer*/) const = 0;

    virtual std::vector<unsigned int> numTmpObjectsRequired() const = 0;
  }; // class LocalOperatorApplication

  template <class LocalMatrixAssemblerType, class MatrixType>
  class LocalMatrixAssemblerApplicationWrapper : public LocalMatrixAssemblerApplication
  {
  public:
    LocalMatrixAssemblerApplicationWrapper(const std::shared_ptr<const LocalMatrixAssemblerType> _localMatrixAssembler,
                                           std::shared_ptr<MatrixType> _matrix)
      : localMatrixAssembler_(_localMatrixAssembler)
      , matrix_(_matrix)
    {
    }

    virtual void apply(const AnsatzSpaceType& _ansatzSpace, const TestSpaceType& _testSpace,
                       const IntersectionType& _intersection, LocalMatricesContainerType& _localMatricesContainer) const
    {
      localMatrixAssembler_->assembleLocal(_intersection, _ansatzSpace, _testSpace, *matrix_, _localMatricesContainer);
    } // virtual void applyLocal(...) const

    virtual std::vector<unsigned int> numTmpObjectsRequired() const
    {
      return localMatrixAssembler_->numTmpObjectsRequired();
    } // virtual std::vector< unsigned int > numTmpObjectsRequired() const

  private:
    const std::shared_ptr<const LocalMatrixAssemblerType> localMatrixAssembler_;
    std::shared_ptr<MatrixType> matrix_;
  }; // class LocalMatrixAssemblerApplicationWrapper

public:
  Boundary(const BoundaryGridPartType& boundaryGridPart, const std::shared_ptr<const BoundaryInfoType> boundaryInfo,
           const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace)
    : boundaryGridPart_(boundaryGridPart)
    , boundaryInfo_(boundaryInfo)
    , ansatzSpace_(ansatzSpace)
    , testSpace_(testSpace)
  {
  }

  ~Boundary()
  {
    for (auto& localMatrixAssembler : localMatrixAssemblers_)
      delete localMatrixAssembler;
    //    for (auto& localVectorAssembler: localVectorAssemblers_)
    //      delete localVectorAssembler;
  }

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

  const std::shared_ptr<const BoundaryInfoType> boundaryInfo() const
  {
    return boundaryInfo_;
  }

  template <class LocalMatrixAssemblerType, class MatrixType>
  void addLocalMatrixAssembler(const std::shared_ptr<const LocalMatrixAssemblerType> _localMatrixAssembler,
                               std::shared_ptr<MatrixType> _matrix)
  {
    typedef LocalMatrixAssemblerApplicationWrapper<LocalMatrixAssemblerType, MatrixType> WrapperType;
    WrapperType* wrapper = new WrapperType(_localMatrixAssembler, _matrix);
    localMatrixAssemblers_.push_back(wrapper);
  }

  //  template< class DirichletAssemblerType, class MatrixBackendType >
  void assemble(/*const DirichletAssemblerType& dirichletAssembler,
                MatrixBackendType& matrix*/) const
  {
    // preparations
    typedef typename BoundaryGridPartType::template Codim<0>::EntityType EntityType;
    typedef typename IntersectionType::EntityPointer EntityPointerType;
    // common tmp storage for all entities
    // * for the matrix assemblers
    std::vector<unsigned int> numberOfTmpMatricesNeeded(2, 0);
    for (unsigned int ii = 0; ii < localMatrixAssemblers_.size(); ++ii) {
      const std::vector<unsigned int> tmp = localMatrixAssemblers_[ii]->numTmpObjectsRequired();
      numberOfTmpMatricesNeeded[0]        = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
      numberOfTmpMatricesNeeded[1]        = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
    }
    std::vector<LocalMatrixType> tmpLocalAssemblerMatrices(
        numberOfTmpMatricesNeeded[0],
        LocalMatrixType(ansatzSpace_.map().maxLocalSize(), testSpace_.map().maxLocalSize(), RangeFieldType(0.0)));
    std::vector<LocalMatrixType> tmpLocalOperatorMatrices(
        numberOfTmpMatricesNeeded[1],
        LocalMatrixType(ansatzSpace_.map().maxLocalSize(), testSpace_.map().maxLocalSize(), RangeFieldType(0.0)));
    std::vector<std::vector<LocalMatrixType>> tmpLocalMatricesContainer;
    tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
    tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
    // walk the boundary grid part
    for (typename BoundaryGridPartType::template Codim<0>::IteratorType entityIt =
             boundaryGridPart_.template begin<0>();
         entityIt != boundaryGridPart_.template end<0>();
         ++entityIt) {
      const EntityType& entity = *entityIt;
      // walk the intersections
      for (IntersectionIteratorType interectionIt = boundaryGridPart_.ibegin(entity);
           interectionIt != boundaryGridPart_.iend(entity);
           ++interectionIt) {
        // get the intersection
        const IntersectionType& intersection = *interectionIt;
        // assemble local matrices
        for (unsigned int ii = 0; ii < localMatrixAssemblers_.size(); ++ii)
          localMatrixAssemblers_[ii]->apply(testSpace_, ansatzSpace_, intersection, tmpLocalMatricesContainer);
      } // walk the intersections
    } // walk the coupling grid part
  } // void assemble(...) const

private:
  const BoundaryGridPartType& boundaryGridPart_;
  const std::shared_ptr<const BoundaryInfoType> boundaryInfo_;
  const AnsatzSpaceType& ansatzSpace_;
  const TestSpaceType& testSpace_;
  std::vector<LocalMatrixAssemblerApplication*> localMatrixAssemblers_;
  //  std::vector< LocalVectorAssemblerApplication* > localVectorAssemblers_;
}; // class Boundary

} // namespace Multiscale

} // namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_BOUNDARY_HH
