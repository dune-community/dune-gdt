#ifndef DUNE_FUNCTIONALS_CONTAINER_SPARSITYPATTERN_HH
#define DUNE_FUNCTIONALS_CONTAINER_SPARSITYPATTERN_HH

#include <vector>
#include <set>

namespace Dune {

namespace Functionals {

namespace Container {
/**
 * @class SparsityPattern
 * Class for storing the sparsity pattern for a sparse matrix
 */
class SparsityPattern
{
public:
  //! type for saving the sparsity pattern
  typedef std::vector<std::set<unsigned int>> SparsityPatternContainerType;

  //! type for iterating through a row
  typedef std::set<unsigned int>::const_iterator NonZeroColIterator;


  /**
   * @brief constructor
   *
   * @param rowSize number of rows for the sparsity pattern
   */
  SparsityPattern(unsigned int rowSize)
    : sparsityPattern_(rowSize)
    , sizeN_(rowSize)
  {
  }

  /**
   * @brief insert a position (row, col) for a nonzero entry
   */
  void insert(unsigned int row, unsigned int col)
  {
    sparsityPattern_[row].insert(col);
  }

  /**
   * @brief remove a position (row, col) for a nonzero entry
   */
  void erase(unsigned int row, unsigned int col)
  {
    sparsityPattern_[row].erase(col);
  }

  /**
   * @brief checks whether block is zero
   */
  bool isZero(unsigned int row, unsigned int col)
  {
    return (sparsityPattern_[row].count(col) == 0);
  }

  /**
   * @brief number of nonzeros in row "row" (counted in blocks)
   */
  unsigned int countNonZeros(unsigned int row)
  {
    return sparsityPattern_[row].size();
  }

  /**
   * @brief number of rows (counted in blocks)
   */
  unsigned int size()
  {
    return sizeN_;
  }

  /**
   * @brief number of rows (counted in blocks)
   */
  unsigned int N()
  {
    return sizeN_;
  }

  /**
   *  @brief get pointer to the first nonzero entry
   *  in row "row"
   */
  NonZeroColIterator begin(unsigned int row)
  {
    sparsityPattern_[row].begin();
  }

  /**
   *  @brief get pointer pointing behind the last nonzero entry
   *  in row "row"
   */
  NonZeroColIterator end(unsigned int row)
  {
    sparsityPattern_[row].end();
  }

private:
  SparsityPatternContainerType sparsityPattern_;
  unsigned int sizeN_;
};

/**
 * @class DefaultSparsityPattern
 * This default implementation sets the diagonal elements
 * to nonzero.
 */
class DefaultSparsityPattern : public SparsityPattern
{

  /**
   * @brief constructor
   */
  DefaultSparsityPattern(unsigned int rowSize)
    : SparsityPattern(rowSize)
  {
    for (int i = 0; i < rowSize; ++i)
      insert(i, i);
  }

}; // end of class DefaultSparsityPattern

} // end of namespace Container

} // end of namespace Functionals

} // end of namespace Dune

#endif // DUNE_FUNCTIONALS_CONTAINER_SPARSITYPATTERN_HH
