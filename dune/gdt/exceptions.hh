// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016)

#ifndef DUNE_GDT_EXCEPTIONS_HH
#define DUNE_GDT_EXCEPTIONS_HH

#include <dune/xt/common/exceptions.hh>

namespace Dune {
namespace GDT {
namespace Exceptions {


class operator_error : public Dune::Exception
{
};

class local_operator_error : public operator_error
{
};

class numerical_flux_error : public local_operator_error
{
};

class prolongation_error : public operator_error
{
};

class projection_error : public operator_error
{
};

class space_error : public Dune::Exception
{
};

class restricted_space_error : public space_error
{
};


} // namespace Exceptions
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EXCEPTIONS_HH
