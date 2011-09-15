#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH

// std includes
#include <vector>

namespace Dune {

namespace DetailedDiscretizations {

namespace Assembler {

namespace Local {

namespace Codim0 {

template <class LocalFunctionalImp>
class Vector
{
public:
  typedef LocalFunctionalImp LocalFunctionalType;

  typedef Vector<LocalFunctionalType> ThisType;

  typedef typename LocalFunctionalType::RangeFieldType RangeFieldType;

  //! constructor
  Vector(const LocalFunctionalType localFunctional)
    : localFunctional_(localFunctional)
  {
  }

  //! copy constructor
  Vector(const ThisType& other)
    : localFunctional_(other.localFunctional())
  {
  }

  const LocalFunctionalType& localFunctional() const
  {
    return localFunctional_;
  }

  template <class TestSpaceType, class EntityType, class VectorType, class LocalVectorType>
  void assembleLocal(const TestSpaceType& testSpace, const EntityType& entity, VectorType& vector,
                     std::vector<LocalVectorType>& tmpLocalVectors) const
  {
    // get the basefunctionset
    typedef typename TestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalTestBaseFunctionSetType;

    const LocalTestBaseFunctionSetType localTestBaseFunctionSet = testSpace.baseFunctionSet().local(entity);

    // check, if we have enough tmp local vectors
    if (tmpLocalVectors.size() < 1) {
      tmpLocalVectors.resize(1, LocalVectorType(testSpace.map().maxLocalSize(), RangeFieldType(0.0)));
    }

    // write local functional application to tmpLocalVector
    localFunctional_.applyLocal(localTestBaseFunctionSet, tmpLocalVectors[0]);

    // write local vector to global
    addToVector(testSpace, entity, tmpLocalVectors[0], vector);
  }

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  template <class TestSpaceType, class EntityType, class LocalVectorType, class VectorType>
  void addToVector(const TestSpaceType& testSpace, const EntityType& entity, const LocalVectorType& localVector,
                   VectorType& vector) const
  {
    for (unsigned int j = 0; j < testSpace.baseFunctionSet().local(entity).size(); ++j) {
      const unsigned int globalJ = testSpace.map().toGlobal(entity, j);

      vector[globalJ] += localVector[j];
    }
  } // end method addToVector

  const LocalFunctionalType localFunctional_;

}; // end class Vector

} // end namespace Codim0

} // end namespace Local

} // end namespace Assembler

} // end namespace DetailedDiscretization

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH
