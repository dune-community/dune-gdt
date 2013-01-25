#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_VECTOR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_VECTOR_HH

#include <vector>

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/common/vector.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace Assembler {
namespace Local {
namespace Codim1 {
namespace Vector {

template <class LocalFunctionalImp, class BoundaryInfoImp>
class Neumann
{
public:
  typedef LocalFunctionalImp LocalFunctionalType;

  typedef BoundaryInfoImp BoundaryInfoType;

  typedef Neumann<LocalFunctionalType, BoundaryInfoType> ThisType;

  Neumann(const LocalFunctionalType& _localFunctional, const Dune::shared_ptr<const BoundaryInfoType> _boundaryInfo)
    : localFunctional_(_localFunctional)
    , boundaryInfo_(_boundaryInfo)
  {
  }

  const LocalFunctionalType& localFunctional() const
  {
    return localFunctional_;
  }

  const Dune::shared_ptr<const BoundaryInfoType> boundaryInfo() const
  {
    return boundaryInfo_;
  }

  std::vector<unsigned int> numTmpObjectsRequired() const
  {
    return {1, localFunctional_.numTmpObjectsRequired()};
  } // std::vector< unsigned int > numTmpObjectsRequired() const

  template <class TestSpaceType, class EntityType, class VectorType, class LocalVectorType>
  void assembleLocal(const TestSpaceType& testSpace, const EntityType& entity, VectorType& vector,
                     std::vector<std::vector<LocalVectorType>>& tmpLocalVectorsContainer) const
  {
    // get the local basefunction set
    typedef typename TestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalTesBaseFunctionSetType;
    const LocalTesBaseFunctionSetType localTestBaseFunctionSet = testSpace.baseFunctionSet().local(entity);

    // check tmp local vectors
    typedef typename LocalFunctionalType::FunctionSpaceType::RangeFieldType RangeFieldType;
    assert(tmpLocalVectorsContainer.size() > 1);
    std::vector<LocalVectorType>& tmpLocalVectors = tmpLocalVectorsContainer[0];
    if (tmpLocalVectors.size() < numTmpObjectsRequired()[0])
      tmpLocalVectors.resize(numTmpObjectsRequired()[0],
                             LocalVectorType(testSpace.map().maxLocalSize(), RangeFieldType(0)));

    // do loop over all intersections
    typedef typename TestSpaceType::GridViewType GridViewType;
    typedef typename GridViewType::IntersectionIterator IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename IntersectionType::EntityPointer EntityPointerType;
    const GridViewType& gridView = testSpace.gridView();
    for (IntersectionIteratorType intersectionIt = gridView.ibegin(entity); intersectionIt != gridView.iend(entity);
         ++intersectionIt) {
      const IntersectionType& intersection = *intersectionIt;
      if (boundaryInfo_->neumann(intersection)) {
        Dune::Stuff::Common::clear(tmpLocalVectors[0]);
        localFunctional_.applyLocal(
            localTestBaseFunctionSet, intersection, tmpLocalVectors[0], tmpLocalVectorsContainer[1]);
        addToVector(testSpace, entity, tmpLocalVectors[0], vector);
      }
    } // do loop over all intersections
  } // void assembleLocal(...)

private:
  Neumann(const ThisType&);
  ThisType& operator=(const ThisType&);

  template <class TestSpaceType, class EntityType, class LocalVectorType, class SystemVectorType>
  void addToVector(const TestSpaceType& testSpace, const EntityType& entity, const LocalVectorType& localVector,
                   SystemVectorType& systemVector) const
  {
    for (unsigned int j = 0; j < testSpace.baseFunctionSet().local(entity).size(); ++j) {
      const unsigned int globalJ = testSpace.map().toGlobal(entity, j);
      systemVector.add(globalJ, localVector[j]);
    }
  } // vodi addToVector(...)

  const LocalFunctionalType& localFunctional_;
  const Dune::shared_ptr<const BoundaryInfoType> boundaryInfo_;
}; // class Neumann

} // namespace Vector
} // namespace Codim1
} // namespace Local
} // namespace Assembler
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_VECTOR_HH
