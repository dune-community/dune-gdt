// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016 - 2018)

/* begin dune-gdt */
// NEVER delete/alter above comment, dune's cmake relies on it

/* Define to the version of dune-gdt */
#define DUNE_GDT_VERSION ${DUNE_GDT_VERSION}

/* Define to the major version of dune-gdt */
#define DUNE_GDT_VERSION_MAJOR ${DUNE_GDT_VERSION_MAJOR}

/* Define to the minor version of dune-gdt */
#define DUNE_GDT_VERSION_MINOR ${DUNE_GDT_VERSION_MINOR}

/* Define to the revision of dune-gdt */
#define DUNE_GDT_VERSION_REVISION ${DUNE_GDT_VERSION_REVISION}

// alberta and lpsolve both define a clashing get_max_level
#if HAVE_LPSOLVE && HAVE_ALBERTA
#undef HAVE_LPSOLVE
#endif

//! TODO make this a cmake cache var
#define GDT_BINDINGS_GRID YASP_2D_EQUIDISTANT_OFFSET
/* end dune-gdt */
// NEVER delete/alter above comment, dune's cmake relies on it
