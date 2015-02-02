// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EXCEPTIONS_HH
#define DUNE_GDT_EXCEPTIONS_HH

#include <dune/stuff/common/exceptions.hh>

namespace Dune {
namespace GDT {
namespace Exceptions {


class operator_error : public Dune::Exception
{
};

class prolongation_error : operator_error
{
};

class darcy_operator_error : public operator_error
{
};


} // namespace Exceptions
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EXCEPTIONS_HH
