// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include "projections/lagrange.hh"
#include "spaces/cg/fem.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_FEM


typedef testing::Types<SPACES_CG_FEM(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(LagrangeProjectionOperatorTest, SpaceTypes);
TYPED_TEST(LagrangeProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LagrangeProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LagrangeProjectionOperatorTest, free_function_callable)
{
  this->free_function_callable();
}
TYPED_TEST(LagrangeProjectionOperatorTest, produces_correct_results)
{
  this->produces_correct_results();
}


#else // HAVE_DUNE_FEM


TEST(DISABLED_LagrangeProjectionOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, free_function_callable)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_FEM
// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include "projections/lagrange.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(LagrangeProjectionOperatorTest, SpaceTypes);
TYPED_TEST(LagrangeProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LagrangeProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LagrangeProjectionOperatorTest, free_function_callable)
{
  this->free_function_callable();
}
TYPED_TEST(LagrangeProjectionOperatorTest, produces_correct_results)
{
  this->produces_correct_results();
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_LagrangeProjectionOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, free_function_callable)
{
}
TEST(DISABLED_LagrangeProjectionOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include "projections/lagrange.hh"
#include "spaces/cg/fem.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_FEM


typedef testing::Types<SPACES_CG_FEM(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(LagrangeProjectionLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, produces_correct_results)
{
  this->produces_correct_results(this->default_tolerance);
}


#else // HAVE_DUNE_FEM


TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_FEM
// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include "projections/lagrange.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(LagrangeProjectionLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(LagrangeProjectionLocalizableOperatorTest, produces_correct_results)
{
  this->produces_correct_results(this->default_tolerance);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_LagrangeProjectionLocalizableOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/type_traits.hh>

#include "projections/l2-global.hh"
#include "spaces/cg/fem.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_FEM


typedef testing::Types<SPACES_CG_FEM(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2GlobalProjectionOperatorTest, SpaceTypes);
TYPED_TEST(L2GlobalProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2GlobalProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2GlobalProjectionOperatorTest, produces_correct_results)
{
  typedef typename TypeParam::GridViewType::Grid Grid;
  const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}


#else // HAVE_DUNE_FEM


TEST(DISABLED_L2GlobalProjectionOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2GlobalProjectionOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2GlobalProjectionOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_FEM
// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx>

#include "projections/l2-global.hh"
#include "spaces/cg/pdelab.hh"
using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2GlobalProjectionOperatorTest, SpaceTypes);
TYPED_TEST(L2GlobalProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2GlobalProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2GlobalProjectionOperatorTest, produces_correct_results)
{
  typedef typename TypeParam::GridViewType::Grid Grid;
  const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2GlobalProjectionOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2GlobalProjectionOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2GlobalProjectionOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/type_traits.hh>

#include "projections/l2-global.hh"
#include "spaces/cg/fem.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_FEM


typedef testing::Types<SPACES_CG_FEM(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2GlobalProjectionLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(L2GlobalProjectionLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2GlobalProjectionLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2GlobalProjectionLocalizableOperatorTest, produces_correct_results)
{
  typedef typename TypeParam::GridViewType::Grid Grid;
  const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}


#else // HAVE_DUNE_FEM


TEST(DISABLED_L2GlobalProjectionLocalizableOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2GlobalProjectionLocalizableOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2GlobalProjectionLocalizableOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_FEM
// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx>

#include "projections/l2-global.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    SpaceTypes;

TYPED_TEST_CASE(L2GlobalProjectionLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(L2GlobalProjectionLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2GlobalProjectionLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2GlobalProjectionLocalizableOperatorTest, produces_correct_results)
{
  typedef typename TypeParam::GridViewType::Grid Grid;
  const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2GlobalProjectionLocalizableOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2GlobalProjectionLocalizableOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2GlobalProjectionLocalizableOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
