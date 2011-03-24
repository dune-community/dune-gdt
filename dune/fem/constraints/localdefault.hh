#ifndef DUNE_FEM_FUNCTIONALS_CONSTRAINTS_LOCALDEFAULT_HH
#define DUNE_FEM_FUNCTIONALS_CONSTRAINTS_LOCALDEFAULT_HH


namespace Dune
{

namespace Functionals
{

namespace Constraints
{

template<typename FieldType, int maxRows, int maxColumns>
class LocalDefault
{
private:
  typedef FieldVector<unsigned int, maxRows>
    RowDofs;

  typedef FieldVector<unsigned int, maxColumns>
    ColumnDofs;

  typedef FieldMatrix<FieldType, maxRows, maxColumns>
    MatrixType;

public:
  LocalDefault( int numColumns=0 )
    : rowDofs_(0),
      columnDofs_(0),
      matrix_(0),
      numRows_(0),
      numColumns_(numColumns)
  {
  }

  void setRowDofsSize( int numRows )
  {
    assert(numRows < maxRows);
    numRows_ = numRows;
  }

  void setColumnDofsSize( int numColumns )
  {
    assert(numColumns < maxColumns);
    numColumns_ = numColumns;
  }

  unsigned int rowDofsSize() const
  {
    return numRows_;
  }

  unsigned int columnDofsSize() const
  {
    return numColumns_;
  }

  void setRowDofs( unsigned int i,
                   unsigned int globalDof )
  {
    rowDofs_[i] = globalDof;
  }

  void setColumnDofs( unsigned int i,
                   unsigned int globalDof )
  {
    columnDofs_[i] = globalDof;
  }

  unsigned int rowDofs( unsigned int i ) const
  {
    assert(i < maxRows);
    return rowDofs_[i];
  }

  unsigned int columnDofs( unsigned int i ) const
  {
    assert(i < maxColumns);
    return columnDofs_[i];
  }

  void setLocalMatrix( unsigned int i,
                       unsigned int j,
                       FieldType val )
  {
    matrix_[i][j] = val;
  }

  FieldType localMatrix( unsigned int i,
                         unsigned int j ) const
  {
    assert(i < maxRows);
    assert(i < maxColumns);
    return matrix_[i][j];
  }

private:
  RowDofs      rowDofs_;
  ColumnDofs   columnDofs_;
  MatrixType   matrix_;
  unsigned int numRows_;
  unsigned int numColumns_;

}; // end class Dirichlet

} // end of namespace Constraints

} // end of namespace Functionals

} // end of namespace Dune

#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_CONSTRAINTS_LOCALDEFAULT_HH */
