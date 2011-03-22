#ifndef DUNE_FUNCTIONALS_SPARSITYPATTERN
#define DUNE_FUNCTIONALS_SPARSITYPATTERN

#include <vector>
#include <set>

namespace Dune
{
namespace Functionals
{

  /**
   * @class SparsityPattern
   * Class for storing the sparsity pattern for a sparse matrix
   */
  class SparsityPattern
  {
  public:
    typedef std::vector< std::set< unsigned int > >
      SparsityPatternContainerType;
   
    /** 
     * @brief constructor
     *  
     * @param rowSize number of rows for the sparsity pattern
     */
    SparsityPattern( unsigned int rowSize )
      : sparsityPattern_( rowSize ), sizeN_( rowSize ) 
    {
    }

    /**
     * @brief insert a position (row, col) for a nonzero entry
     */
    void insert( unsigned int row, unsigned int col )
    {
      sparsityPattern_[row].insert( col );
    }

    /**
     * @brief remove a position (row, col) for a nonzero entry
     */
    void erase( unsigned int row, unsigned int col )
    {
      sparsityPattern_[row].erase( col );
    }

    /**
     * @brief checks whether block is zero
     */
    bool isZero( unsigned int row, unsigned int col )
    {
      return ( sparsityPattern_[row].count( col ) == 0 )
    }

    /**
     * @brief number of nonzeros in row "row" (counted in blocks)
     */
    unsigned int countNonZeros( unsigned int row )
    {
      return sparsityPattern_[row].size();
    }

    /**
     * @brief number of rows (counted in blocks)
     */
    unsigned int N()
    { 
      return sizeN_;
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
  class DefaultSparsityPattern
    : public SparsityPattern
  {
     
    /**
     * @brief constructor 
     */
    DefaultSparsityPattern( unsigned int rowSize )
      : SparsityPattern( rowSize )
    {
      for( int i = 0; i < rowSize; ++i )
        insert( i, i );
    }

  }; //end of class DefaultSparsityPattern

} // end of namespace Functionals
} // end of namespace Dune

#endif// DUNE_FUNCTIONALS_SPARSITYPATTERN

