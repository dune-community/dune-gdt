# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2017)

import itertools

from dune.xt.codegen import typeid_to_typedef_name
from grids import LeafGrids, LevelGrids


def CG(cache, base=LeafGrids, orders=range(1,2)):
    cg = base(cache)
    cg.grids = [cg.yasp_view_fmt.format(1)] + [s.format(d) for s, d in itertools.product(cg.all_views_fmt, cg.world_dim)]
    cg.spaces = ['Dune::GDT::ContinuousLagrangeSpace<{}, {}>'.format(grid, order)
              for grid, order in itertools.product(cg.grids, orders)]
    cg.names = [typeid_to_typedef_name(sp) for sp in cg.spaces]
    return cg


def DG(cache, base=LeafGrids):
    dg = base(cache)
    dg.grids = [dg.yasp_view_fmt.format(1)] + [s.format(d) for s, d in itertools.product(dg.all_views_fmt, dg.world_dim)]
    dg.spaces = ['Dune::GDT::DiscontinuousLagrangeSpace<{}, 1, double>'.format(grid)
                 for grid in dg.grids]
    dg.names = [typeid_to_typedef_name(sp) for sp in dg.spaces]
    return dg


def FV(cache, base=LeafGrids, rdim=None):
    fv = base(cache)
    cg = CG(cache, base=base)
    fv.rdim = rdim or [1, 2, 3]
    fv.grids = cg.grids
    fv.spaces = ['Dune::GDT::FvSpace<{}, double, {}>'.format(g, d) for g, d in itertools.product(fv.grids, fv.rdim)]
    fv.names = [typeid_to_typedef_name(sp) for sp in fv.spaces]
    return fv


def RT(cache, base=LeafGrids):
    rt = base(cache)
    rt.spaces = ['Dune::GDT::DunePdelabRtSpaceWrapper<{}, 0, double, {}>'.format(s.format(d), d)
                   for s, d in itertools.product(rt.all_views_fmt, rt.world_dim)]
    rt.names = [typeid_to_typedef_name(sp) for sp in rt.spaces]
    return rt


if __name__ == '__dxt_codegen__':
    # this is executed from spaces.tpl itself
    cg = CG(cache)
    dg = DG(cache)
    fv = FV(cache)
    rt = RT(cache)
    spaces = cg.spaces + dg.spaces + fv.spaces + rt.spaces
    names = cg.names + dg.names + fv.names + rt.names
    spaces_with_names = zip(spaces, names)
