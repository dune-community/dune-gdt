// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH
#define DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH

#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace Dune {
namespace GDT {
namespace TmpStorageProvider {


template< class FieldType = double >
class Matrices
{
  static size_t safely_get(const std::vector< size_t >& vec, const size_t ii)
  {
    assert(ii < vec.size());
    return vec[ii];
  } // ... safely_get(...)

public:
  typedef DynamicMatrix< FieldType > LocalMatrixType;

  Matrices(const std::vector< size_t >& num_tmp_objects,
           const size_t max_rows,
           const size_t max_cols)
    : matrices_({ std::vector< LocalMatrixType >(safely_get(num_tmp_objects, 0),
                                                 LocalMatrixType(max_rows, max_cols, FieldType(0)))
                , std::vector< LocalMatrixType >(safely_get(num_tmp_objects, 1),
                                                 LocalMatrixType(max_rows, max_cols, FieldType(0))) })
    , indices_(4, Dune::DynamicVector< size_t >(std::max(max_rows, max_cols)))
  {}

  virtual ~Matrices() {}

  std::vector< std::vector< LocalMatrixType > >& matrices()
  {
    return matrices_;
  }

  std::vector< Dune::DynamicVector< size_t > >& indices()
  {
    return indices_;
  }

protected:
  std::vector< std::vector< LocalMatrixType > > matrices_;
  std::vector< DynamicVector< size_t > > indices_;
}; // class Matrices


template< class FieldType = double >
class Vectors
{
  static size_t safely_get(const std::vector< size_t >& vec, const size_t ii)
  {
    assert(ii < vec.size());
    return vec[ii];
  } // ... safely_get(...)

public:
  typedef DynamicVector< FieldType > LocalVectorType;

  Vectors(const std::vector< size_t >& num_tmp_objects,
          const size_t max_size)
    : vectors_({ std::vector< LocalVectorType >(safely_get(num_tmp_objects, 0),
                                                LocalVectorType(max_size, FieldType(0)))
               , std::vector< LocalVectorType >(safely_get(num_tmp_objects, 1),
                                                LocalVectorType(max_size, FieldType(0))) })
    , indices_(max_size)
  {}

  virtual ~Vectors() {}

  std::vector< std::vector< LocalVectorType > >& vectors()
  {
    return vectors_;
  }

  Dune::DynamicVector< size_t >& indices()
  {
    return indices_;
  }

protected:
  std::vector< std::vector< LocalVectorType > > vectors_;
  Dune::DynamicVector< size_t > indices_;
}; // class Vectors


} // namespace TmpStorageProvider
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH
