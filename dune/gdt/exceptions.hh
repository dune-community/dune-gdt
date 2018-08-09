// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016)

#ifndef DUNE_GDT_EXCEPTIONS_HH
#define DUNE_GDT_EXCEPTIONS_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/grid/exceptions.hh>
#include <dune/xt/functions/exceptions.hh>

namespace Dune {
namespace GDT {
namespace Exceptions {


class not_bound_to_an_element_yet : public XT::Grid::Exceptions::not_bound_to_an_element_yet
{
};

class dof_vector_error : public Exception
{
};

class integrand_error : public Exception
{
};

class assembler_error : public Exception
{
};

class functional_error : public Exception
{
};

class numerical_flux_error : public Exception
{
};

class operator_error : public Exception
{
};

class prolongation_error : public operator_error
{
};

class projection_error : public operator_error
{
};

class finite_element_error : public Exception
{
};

class space_error : public Exception
{
};

class restricted_space_error : public space_error
{
};

class mapper_error : public space_error
{
};

class basis_error : public space_error
{
};

class discrete_function_error : public Exception
{
};


} // namespace Exceptions
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EXCEPTIONS_HH
