/**
  \file matrix.hh
  \todo doc
  **/

#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_COMBINED_MATRIX_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_COMBINED_MATRIX_HH

// std includes
#include <vector>

namespace Dune {

namespace DetailedDiscretizations {

namespace Assembler {

namespace Local {

namespace Combined {

/**
  \brief  Combined local matrix assembler.
  \todo   doc
  **/
template <class FirstImp, class SecondImp>
class Matrix
{
public:
  typedef FirstImp FirstType;

  typedef SecondImp SecondType;

  typedef Matrix<FirstType, SecondType> ThisType;

  /**
    \brief  Constructor.
    \todo   doc
    **/
  Matrix(const FirstType& first, const SecondType& second)
    : first_(first)
    , second_(second)
  {
  }

private:
  /**
    \brief  Copy constructor.
    \todo   doc
    **/
  Matrix(const ThisType& other)
    : first_(other.first())
    , second_(other.second())
  {
  }

public:
  /**
    \brief  Returns first local matrix assembler.
    \todo   doc
    **/
  const FirstType& first() const
  {
    return first_;
  }

  /**
    \brief  Returns second local matrix assembler.
    \todo   doc
    **/
  const SecondType& second() const
  {
    return second_;
  }

  /**
    \brief  Returns the number of max tmp storage, which are required.
    \todo   doc
    **/
  std::vector<unsigned int> numTmpObjectsRequired() const
  {
    std::vector<unsigned int> ret(2, 0);
    std::vector<unsigned int> numTmpObjectsRequiredFirst  = first_.numTmpObjectsRequired();
    std::vector<unsigned int> numTmpObjectsRequiredSecond = second_.numTmpObjectsRequired();
    ret[0]                                                = std::max(numTmpObjectsRequiredFirst[0], numTmpObjectsRequiredSecond[0]);
    ret[1]                                                = std::max(numTmpObjectsRequiredFirst[1], numTmpObjectsRequiredSecond[1]);
    return ret;
  }

  /**
    \brief  Main method.
    \todo   doc
    **/
  template <class AnsatzSpaceType, class TestSpaceType, class EntityType, class SystemMatrixType, class LocalMatrixType>
  void assembleLocal(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace, const EntityType& entity,
                     SystemMatrixType& systemMatrix,
                     std::vector<std::vector<LocalMatrixType>>& tmpLocalMatricesContainer) const
  {
    first_.assembleLocal(ansatzSpace, testSpace, entity, systemMatrix, tmpLocalMatricesContainer);
    second_.assembleLocal(ansatzSpace, testSpace, entity, systemMatrix, tmpLocalMatricesContainer);
  } // end method assembleLocal

private:
  /**
    \brief  Assignment operator.
    \todo   doc
    **/
  ThisType& operator=(const ThisType&);

  const FirstType& first_;
  const SecondType& second_;
}; // end class Matrix

} // end namespace Combined

} // end namespace Local

} // end namespace Assembler

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_COMBINED_MATRIX_HH
