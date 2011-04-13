#ifndef DUNE_FEM_FUNCTIONALS_COMMON_LOCALMATRIX_HH
#define DUNE_FEM_FUNCTIONALS_COMMON_LOCALMATRIX_HH

// dunefem-functionals include
#include "localvector.hh"

namespace Dune {

namespace Functionals {

namespace Common {

/**
  * \brief  This class represents a local DoF Matrix.
  *         It is based upon std::vector and should be replaced by something clever in the future!
  *
  * \tparam ElementType
  *         Type of one element, usually double or RangeFieldType.
  *
  * \todo   Doc me, please!
  **/
template <class ElementType>
class LocalMatrix
{
public:
  typedef LocalVector<ElementType> RowType;

  /**
    * \brief      Initializes an empty matrix, according to the given size.
    *
    * \param[in]  rows
    *             Number of rows, the matrix will have.
    *
    * \param[in]  cols
    *             Number of collumns, the matrix will have.
    **/
  LocalMatrix(const unsigned int rows, const unsigned int cols)
  {
    // resize
    storage_.resize(rows, RowType(cols));
  }

  /**
    * \brief    Returns the storage size (n*m).
    */
  const unsigned int size() const
  {
    return storage_.size() * storage_[0].size();
  }

  RowType& operator[](const unsigned int i)
  {
    return storage_[i];
  }

  const RowType& operator[](const unsigned int i) const
  {
    return storage_[i];
  }

  /**
    * \brief    Matrix product with a local DoF vector.
    **/
  LocalVector<ElementType> operator*(const LocalVector<ElementType>& other) const
  {
    assert(cols() == other.size());

    LocalVector<ElementType> result(rows());

    for (unsigned int i = 0; i < rows(); ++i) {
      for (unsigned int j = 0; j < cols(); ++j) {
        result[i] += storage_[i][j] * other[j];
      }
    }

    return result;
  }

  /**
   * @brief Return number of rows.
   */
  unsigned int N() const
  {
    return storage_.size();
  }

  /**
   * @brief Return number of rows.
   */
  unsigned int rows() const
  {
    return storage_.size();
  }

  /**
   * @brief Return number of cols.
   */
  unsigned int M() const
  {
    return storage_[0].size();
  }

  /**
   * @brief Return number of cols.
   */
  unsigned int cols() const
  {
    return storage_[0].size();
  }

private:
  std::vector<RowType> storage_;

}; // end class LocalDoFVector

} // end namespace Functionals

} // end namespace Common

} // end namespace Dune

#endif // end of include guard: DUNE_FEM_FUNCTIONALS_COMMON_LOCALMATRIX_HH
