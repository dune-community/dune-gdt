#ifndef DUNE_DETAILED_DISCRETIZATIONS_CONSTRAINTS_LOCALDEFAULT_HH
#define DUNE_DETAILED_DISCRETIZATIONS_CONSTRAINTS_LOCALDEFAULT_HH

// dune-common
#include <dune/common/fmatrix.hh>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

/**
 * @addtogroup Constraints Introduction to Constraints
 *
 * Let @f$V_h@f$ be a linear discrete function space, @f$M \in \mathbb{N}_{>0}@f$.
 * Then a set of linear functionals @f$C=\{ C_1, ..., C_M \}@f$
 * on @f$V_h@f$ with the condition
 * @f{align*}{
 *   C_i[v] = 0 \quad \forall 1 \leq i \leq M
 * @f}
 * is called a @e constraint for @f$v@f$.
 * Thus each linear functional implies a constraint.
 *
 * Let @f$V_h|_{G}@f$ be a restriction to a local domain @f$G@f$,
 * for example a support of a basis function.
 * Then we call
 * @f{align*}{
 *   C|_G = \{C_1|_G,\ldots,C_M|_G\}
 * @f}
 * a @link Constraints::LocalDefault @e local @e constraint @endlink.
 *
 */

/**
 * @brief Contains various constraints.
 *
 * @ingroup Constraints
 */
namespace Constraints
{

/**
 * @brief Class implementing the local constraint interface.
 *
 * The most important method @link LocalDefault::localMatrix()
 * localMatrix() @endlink returns a local matrix. By means of the
 * methods @link LocalDefault::rowDofs() rowDofs() @endlink,
 * @link LocalDefault::columnDofs() columnDofs() @endlink
 * it is possible to get for each row and column number in the
 * local matrix the corresponding row and column number
 * in the global matrix.
 *
 * @tparam FieldType The type of the entries in the local matrix, i.e. double.
 * @tparam maxRows The maximum number of rows of the local matrix.
 * @tparam maxColumns The maximum number of columns of the local matrix.
 *
 * @todo Think about the name of this class!
 *
 * @ingroup Constraints
 */
template< typename FieldType, int maxRows, int maxColumns >
class LocalDefault
{
private:
  static const unsigned int maxRows_ = maxRows;

  static const unsigned int maxCols_ = maxColumns;

  typedef FieldVector< unsigned int, maxRows_ >
    RowDofs;

  typedef FieldVector< unsigned int, maxCols_ >
    ColumnDofs;

  typedef FieldMatrix< FieldType, maxRows_, maxCols_ >
    MatrixType;

public:
  typedef LocalDefault< FieldType, maxRows_, maxCols_ >
    ThisType;

  /**
   * @brief Constructor.
   *
   * @param numColumns The size of the columns for the degrees of
   *        freedom in the local matrix.
   */
  LocalDefault( int numColumns = 0 )
    : rowDofs_( 0.0 ),
      columnDofs_( 0.0 ),
      matrix_( 0.0 ),
      numRows_( 0 ),
      numColumns_( numColumns )
  {
  }

  //! copy constructor
  LocalDefault( const ThisType& other )
    : rowDofs_( 0.0 ),
      columnDofs_( 0.0 ),
      matrix_( 0.0 ),
      numRows_( other.getRowDofsSize() ),
      numColumns_( other.getColumnDofsSize() )
  {
    std::cout << "Constraints::LocalDefault::LocalDefault( const ThisType& )" << std::endl;
    for( unsigned int i = 0; i < rowDofs_.N(); ++i )
    {
      rowDofs_[i] = other.rowDofs( i );
    }
    for( unsigned int i = 0; i < columnDofs_.N(); ++i )
    {
      columnDofs_[i] = other.columnDofs( i );
    }
    for( unsigned int i = 0; i < matrix_.N(); ++i )
    {
      for( unsigned int j = 0; j < matrix_.M(); ++j )
      {
        setLocalMatrix( i, j, other.localMatrix( i, j ) );
      }
    }
  }

  /**
   * @brief Sets the number of rows for the local matrix.
   *
   * @param numRows The number of rows.
   */
  void setRowDofsSize( unsigned int numRows )
  {
    assert( numRows <= maxRows_ );
    numRows_ = numRows;
  }

  unsigned int getRowDofsSize() const
  {
    return numRows_;
  }

  /**
   * @brief Sets the number of columns for the local matrix.
   *
   * @param numRows The number of columns.
   */
  void setColumnDofsSize( unsigned int numColumns )
  {
    assert( numColumns <= maxCols_ );
    numColumns_ = numColumns;
  }

  unsigned int getColumnDofsSize() const
  {
    return numColumns_;
  }

  /**
   * @brief Gets the number of rows for the local matrix.
   *
   * @result The number of rows.
   */
  unsigned int rowDofsSize() const
  {
    return numRows_;
  }

  /**
   * @brief Gets the number of columns for the local matrix.
   *
   * @result The number of columns.
   */
  unsigned int columnDofsSize() const
  {
    return numColumns_;
  }

  /**
   * @brief Saves the row mapping between local matrix
   * and global matrix.
   *
   * @param i The row number in the local matrix.
   * @param globalDof The row number in global matrix.
   */
  void setRowDofs( unsigned int i, unsigned int globalDof )
  {
    assert( i < maxRows_ );
    rowDofs_[i] = globalDof;
  }

  /**
   * @brief Saves the column mapping between local matrix
   * and global matrix.
   *
   * @param i The column number in the local matrix.
   * @param globalDof The column number in global matrix.
   */
  void setColumnDofs( unsigned int i, unsigned int globalDof )
  {
    columnDofs_[i] = globalDof;
  }

  /**
   * @brief Gets the row mapping between local matrix
   * and global matrix.
   *
   * @param i The row number in the local matrix.
   */
  unsigned int rowDofs( unsigned int i ) const
  {
    assert( i < maxRows_ );
    return rowDofs_[i];
  }

  /**
   * @brief Gets the column mapping between local matrix
   * and global matrix.
   *
   * @param i The column number in the local matrix.
   */
  unsigned int columnDofs( unsigned int i ) const
  {
    assert( i < maxCols_ );
    return columnDofs_[i];
  }

  /**
   * @brief Sets an entry in the local matrix.
   *
   * @param i Row number in the local matrix.
   * @param j Column number in the local matrix.
   * @param val Value of the entry (i,j).
   */
  void setLocalMatrix( unsigned int i,
                       unsigned int j,
                       FieldType val )
  {
    matrix_[i][j] = val;
  }

  /**
   * @brief Gets an entry in the local matrix.
   *
   * @param i Row number in the local matrix.
   * @param j Column number in the local matrix.
   * @return Value of the entry (i,j).
   */
  FieldType localMatrix( unsigned int i,
                         unsigned int j ) const
  {
    assert(i < maxRows);
    assert(j < maxColumns);
    return matrix_[i][j];
  }

private:
  //! assignment operator
  ThisType& operator=( const ThisType& );

  RowDofs      rowDofs_;
  ColumnDofs   columnDofs_;
  MatrixType   matrix_;
  unsigned int numRows_;
  unsigned int numColumns_;

}; // end class Dirichlet

} // end of namespace Constraints

} // namespace Discretizations

} // end of namespace Detailed

} // end of namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_CONSTRAINTS_LOCALDEFAULT_HH
