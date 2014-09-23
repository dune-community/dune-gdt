// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH
#define DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH

#include <dune/stuff/common/tmp-storage.hh>

namespace Dune {
namespace GDT {
namespace TmpStorageProvider {

template <class T>
using Matrices = DSC::TmpMatricesStorage<T>;
template <class T>
using Vectors = DSC::TmpVectorsStorage<T>;

} // namespace TmpStorageProvider
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_TMP_STORAGE_HH
