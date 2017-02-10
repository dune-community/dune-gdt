import itertools

from dune.xt.codegen import typeid_to_typedef_name
from grids import LeafGrids


def CG(cache, base=LeafGrids):
    cg = base(cache)
    cg.fem_grids = [cg.yasp_part_fmt.format(1)] + [s.format(d) for s, d in itertools.product(cg.all_parts_fmt, cg.world_dim)]
    cg.pdelab_grids = [cg.yasp_view_fmt.format(1)] + [s.format(d) for s, d in itertools.product(cg.all_views_fmt, cg.world_dim)]
    cg.fem = ['Dune::GDT::DuneFemCgSpaceWrapper<{}, 1, double, 1>'.format(grid) for grid in cg.fem_grids]
    cg.pdelab = ['Dune::GDT::DunePdelabCgSpaceWrapper<{}, 1, double, 1>'.format(grid) for grid in cg.pdelab_grids]
    cg.spaces = cg.fem + cg.pdelab
    cg.names = [typeid_to_typedef_name(sp) for sp in cg.spaces]
    return cg


def DG(cache, base=LeafGrids):
    dg = base(cache)
    cg = CG(cache, base=base)
    dg.fem_grids = cg.fem_grids
    dg.pdelab_grids = [dg.yasp_view_fmt.format(1)] + [s.format(d) for s, d in itertools.product([dg.alu_cube_view_fmt, dg.yasp_view_fmt], dg.world_dim)]
    dg.fem = ['Dune::GDT::DuneFemDgSpaceWrapper<{}, 1, double, 1>'.format(grid) for grid in dg.fem_grids]
    dg.pdelab = ['Dune::GDT::DunePdelabDgSpaceWrapper<{}, 1, double, 1>'.format(grid) for grid in dg.pdelab_grids]
    dg.spaces = dg.fem + dg.pdelab
    dg.names = [typeid_to_typedef_name(sp) for sp in dg.spaces]
    return dg


def FV(cache, base=LeafGrids, rdim=None):
    fv = base(cache)
    cg = CG(cache, base=base)
    fv.rdim = rdim or [1, 2, 3]
    fv.grids = cg.pdelab_grids
    fv.spaces = ['Dune::GDT::FvSpace<{}, double, {}>'.format(g, d) for g, d in itertools.product(fv.grids, fv.rdim)]
    fv.names = [typeid_to_typedef_name(sp) for sp in fv.spaces]
    return fv


def RT(cache, base=LeafGrids):
    rt = base(cache)
    rt.spaces = ['Dune::GDT::DunePdelabRtSpaceWrapper<{}, 0, double, {}>'.format(s.format(d), d)
                   for s, d in itertools.product(rt.all_views_fmt, rt.world_dim)]
    rt.names = [typeid_to_typedef_name(sp) for sp in rt.spaces]
    return rt