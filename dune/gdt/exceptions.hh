// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)

#ifndef DUNE_GDT_EXCEPTIONS_HH
#define DUNE_GDT_EXCEPTIONS_HH

#include <dune/stuff/common/exceptions.hh>

namespace Dune {
namespace GDT {


class operator_error : public Dune::Exception
{
};

class prolongation_error : public operator_error
{
};

class projection_error : public operator_error
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EXCEPTIONS_HH
