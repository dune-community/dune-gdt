#ifndef DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH
#define DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH

// dune-functionals includes
#include <dune/functionals/common/localvector.hh>

namespace Dune {

namespace Functionals {

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
  Vector(const LocalFunctionalType& localFunctional)
    : localFunctional_(localFunctional)
  {
  }

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
                     LocalVectorType& tmpLocalVector) const
  {
    // write local operator application to tmpLocalMatrix
    localFunctional_.applyLocal(testSpace.localBaseFunctionSet(entity), tmpLocalVector);

    // write local matrix to global
    addToVector(testSpace, entity, tmpLocalVector, vector);
  }

private:
  template <class TestSpaceType, class EntityType, class LocalVectorType, class VectorType>
  void addToVector(const TestSpaceType& testSpace, const EntityType& entity, const LocalVectorType& localVector,
                   VectorType& vector) const
  {
    for (unsigned int j = 0; j < testSpace.localBaseFunctionSet(entity).size(); ++j) {
      const unsigned int globalJ = testSpace.map().toGlobal(entity, j);

      vector[globalJ] += localVector[j];
    }
  } // end method addToVector

  const LocalFunctionalType& localFunctional_;

}; // end class Vector

} // end namespace Codim0

} // end namespace Local

} // end namespace Assembler

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH
