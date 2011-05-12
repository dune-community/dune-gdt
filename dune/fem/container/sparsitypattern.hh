#ifndef DUNE_FUNCTIONALS_CONTAINER_SPARSITYPATTERN_HH
#define DUNE_FUNCTIONALS_CONTAINER_SPARSITYPATTERN_HH

#include <vector>
#include <set>

namespace Dune
{

namespace Functionals
{

namespace Container
{
/**
 * @brief Class for storing the sparsity pattern of a sparse matrix.
 *
 * A sparsity pattern stores nonzero entries (or entries which might be nonzero)
 * in a matrix @f$A=(a_{ij})_{i,j} @f$, 
 * i.e. all entries @f$(i,j)@f$ in the matrix with 
 * @f[(a_{ij})_{i,j}\neq 0\quad\forall i\in\{1,\ldots,m\},j\in\{1,\ldots,n\}@f] 
 * where @f$m@f$ is the number of rows and  @f$n@f$ is the number of columns
 * of the matrix @f$A@f$.
 *
 * Normally, we want to use this class for storing overlapping degrees of freedom
 * of local basis function.
 */
class SparsityPattern
{
public:
  //! Type for saving the sparsity pattern.
  typedef std::vector< std::set< unsigned int > >
    SparsityPatternContainerType;

  //! Type for iterating through a row.
  typedef std::set< unsigned int >::const_iterator
    NonZeroColIterator;



  /**
   * @brief Constructor storing the row size.
   *
   * @param rowSize Number of rows for the sparsity pattern.
   */
  SparsityPattern( unsigned int rowSize )
    : sparsityPattern_( rowSize ), sizeN_( rowSize )
  {
  }

  /**
   * @brief Inserts a nonzero entry.
   *
   * @param row The row number for the nonzero entry.
   * @param col The column number for the nonzero entry.
   */
  void insert( unsigned int row, unsigned int col )
  {
    sparsityPattern_[row].insert( col );
  }

  /**
   * @brief Removes a nonzero entry.
   *
   * @param row The row number for the nonzero entry.
   * @param col The column number for the nonzero entry.
   */
  void erase( unsigned int row, unsigned int col )
  {
    sparsityPattern_[row].erase( col );
  }

  /**
   * @brief Checks whether block is zero.
   *
   * @param row The row number.
   * @param col The column number.
   */
  bool isZero( unsigned int row, unsigned int col )
  {
    return ( sparsityPattern_[row].count( col ) == 0 );
  }

  /**
   * @brief Counts the number of nonzeros in a row (counted in blocks).
   *
   * @param row The row number in the matrix, where we want to count the nonzero entries.
   */
  unsigned int countNonZeros( unsigned int row )
  {
    return sparsityPattern_[row].size();
  }

  /**
   * @brief Returns the number of rows (counted in blocks).
   */
  unsigned int size()
  {
    return sizeN_;
  }

  /**
   * @brief Returns the number of rows (counted in blocks).
   */
  unsigned int N()
  {
    return sizeN_;
  }

  /**
   *  @brief Gets pointer to the first nonzero entry .
   *
   *  @param row The row number.
   */
  NonZeroColIterator begin( unsigned int row )
  {
    return sparsityPattern_[row].begin();
  }

  /**
   *  @brief Gets pointer pointing behind the last nonzero entry.
   *  
   *  @param row The row number.
   */
  NonZeroColIterator end( unsigned int row )
  {
    return sparsityPattern_[row].end();
  }

private:
  SparsityPatternContainerType sparsityPattern_;
  unsigned int sizeN_;
};

/**
 * @brief This default implementation sets the diagonal elements
 * to nonzero.
 */
class DefaultSparsityPattern
  : public SparsityPattern
{

  /**
   * @brief Constructor setting the diagonal elements to nonzero.
   *
   * @param rowSize Number of rows for the sparsity pattern.
   */
  DefaultSparsityPattern( unsigned int rowSize )
    : SparsityPattern( rowSize )
  {
    for( unsigned int i = 0; i < rowSize; ++i )
      insert( i, i );
  }

}; //end of class DefaultSparsityPattern

} // end of namespace Container

} // end of namespace Functionals

} // end of namespace Dune

#endif // DUNE_FUNCTIONALS_CONTAINER_SPARSITYPATTERN_HH

