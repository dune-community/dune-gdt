#ifndef DUNE_FEM_FUNCTIONALS_COMMON_LOCALMATRIX_HH
#define DUNE_FEM_FUNCTIONALS_COMMON_LOCALMATRIX_HH

// dunefem-functionals include
#include "localvector.hh"

namespace Dune {

namespace Functionals {

namespace Common {

/**
  * \brief      This class represents a local DoF Matrix.
  *
  *             It is based upon std::vector and should be replaced by something clever in the future!
  *
  * \todo       Doc me, please!
  **/
template <class ElementType>
class LocalMatrix
{
public:
  /**
    * \brief    Initializes an empty matrix, according to the given size.
    **/
  LocalMatrix(const unsigned int rows, const unsigned int cols)
    : rows_(rows)
    , cols_(cols)
  {
    // resize
    storage_.resize(rows_ * cols_, 0.0);
  }

  /**
   * @brief     Resizes the storage size (rows*cols).
   */
  void resize(const unsigned int rows, const unsigned int cols)
  {
    storage_.resize(rows * cols, 0.0);
  }

  /**
    * \brief    Returns the storage size (n*m).
    */
  const unsigned int size() const
  {
    return rows_ * cols_;
  }

  ElementType get(const unsigned int i, const unsigned int j) const
  {
    return storage_[i * cols_ + j];
  }

  void set(const unsigned int i, const unsigned int j, ElementType val)
  {
    storage_[i * cols_ + j] = val;
  }

  /**
    * \brief    Matrix product with a local DoF vector.
    **/
  LocalVector<ElementType> operator*(const LocalVector<ElementType>& other) const
  {
    assert(cols_ == other.size());

    LocalVector<ElementType> result(rows_);

    for (unsigned ii = 0; ii < rows_; ++ii) {

      for (unsigned jj = 0; jj < cols_; ++jj) {

        result[ii] += get(ii, jj) * other[jj];
      }
    }

    return result;
  }

  /**
   * @brief Return number of rows.
   */
  unsigned int N() const
  {
    return rows_;
  }

  /**
   * @brief Return number of cols.
   */
  unsigned int M() const
  {
    return cols_;
  }

private:
  std::vector<ElementType> storage_;
  unsigned int rows_;
  unsigned int cols_;

}; // end class LocalDoFVector

} // end namespace Functionals

} // end namespace Common

} // end namespace Dune

#endif // end of include guard: DUNE_FEM_FUNCTIONALS_COMMON_LOCALMATRIX_HH
