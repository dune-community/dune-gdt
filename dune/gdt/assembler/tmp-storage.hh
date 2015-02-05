// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH
#define DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH

#include <dune/common/deprecated.hh>

#include <dune/stuff/common/tmp-storage.hh>
#include <dune/stuff/aliases.hh>

namespace Dune {
namespace GDT {
namespace TmpStorageProvider {


template< class T >
class
  DUNE_DEPRECATED_MSG("Use Dune::Stuff::Common::TmpMatricesStorage instead (04.02.2015)!")
      Matrices
  : public DSC::TmpMatricesStorage< T >
{
  typedef DSC::TmpMatricesStorage< T > BaseType;

public:
  template< class... Args >
  Matrices(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
}; // class Matrices


template< class T >
class
  DUNE_DEPRECATED_MSG("Use Dune::Stuff::Common::TmpVectorsStorage instead (04.02.2015)!")
      Vectors
  : public DSC::TmpVectorsStorage< T >
{
  typedef DSC::TmpVectorsStorage< T > BaseType;

public:
  template< class... Args >
  Vectors(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
}; // class Vectors


} // namespace TmpStorageProvider
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH
