import itertools

from dune.xt.codegen import typeid_to_typedef_name


class Grids:
    world_dim = [2, 3]

    alu_cube_part_fmt = 'AluCube{}dLeafGridPartType'
    alu_conf_part_fmt = 'AluConform{}dLeafGridPartType'
    yasp_part_fmt = 'Yasp{}dLeafGridPartType'

    alu_cube_view_fmt = 'AluCube{}dLeafGridViewType'
    alu_conf_view_fmt = 'AluConform{}dLeafGridViewType'
    yasp_view_fmt = 'Yasp{}dLeafGridViewType'

    def __init__(self, cache):
        self.all_parts_fmt = [self.alu_conf_part_fmt, self.alu_cube_part_fmt, self.yasp_part_fmt] if cache['dune-alugrid'] else [self.yasp_part_fmt]
        self.all_views_fmt = [self.alu_conf_view_fmt, self.alu_cube_view_fmt, self.yasp_view_fmt] if cache['dune-alugrid'] else [self.yasp_view_fmt]


class CG(Grids):
    def __init__(self, cache):
        super(CG, self).__init__(cache)
        self.fem_grids = [self.yasp_part_fmt.format(1)] + [s.format(d) for s, d in itertools.product(self.all_parts_fmt, self.world_dim)]
        self.pdelab_grids = [self.yasp_view_fmt.format(1)] + [s.format(d) for s, d in itertools.product(self.all_views_fmt, self.world_dim)]
        self.fem = ['Dune::GDT::DuneFemCgSpaceWrapper<{}, 1, double, 1>'.format(grid) for grid in self.fem_grids]
        self.pdelab = ['Dune::GDT::DunePdelabCgSpaceWrapper<{}, 1, double, 1>'.format(grid) for grid in self.pdelab_grids]
        self.spaces = self.fem + self.pdelab
        self.names = [typeid_to_typedef_name(sp) for sp in self.spaces]


class DG(Grids):
    def __init__(self, cache):
        super(DG, self).__init__(cache)
        cg = CG(cache)
        self.fem_grids = cg.fem_grids
        self.pdelab_grids = [self.yasp_view_fmt.format(1)] + [s.format(d) for s, d in itertools.product([self.alu_cube_view_fmt, self.yasp_view_fmt], self.world_dim)]
        self.fem = ['Dune::GDT::DuneFemDgSpaceWrapper<{}, 1, double, 1>'.format(grid) for grid in self.fem_grids]
        self.pdelab = ['Dune::GDT::DunePdelabDgSpaceWrapper<{}, 1, double, 1>'.format(grid) for grid in self.pdelab_grids]
        self.spaces = self.fem + self.pdelab
        self.names = [typeid_to_typedef_name(sp) for sp in self.spaces]


class FV(Grids):
    def __init__(self, cache):
        super(FV, self).__init__(cache)
        cg = CG(cache)
        self.rdim = [1, 2, 3]
        self.grids = cg.pdelab_grids
        self.spaces = ['Dune::GDT::FvSpace<{}, double, {}>'.format(g, d) for g, d in itertools.product(self.grids, self.rdim)]
        self.names = [typeid_to_typedef_name(sp) for sp in self.spaces]


class RT(Grids):
    def __init__(self, cache):
        super(RT, self).__init__(cache)
        self.spaces = ['Dune::GDT::DunePdelabRtSpaceWrapper<{}, 0, double, {}>'.format(s.format(d), d)
                       for s, d in itertools.product(self.all_views_fmt, self.world_dim)]
        self.names = [typeid_to_typedef_name(sp) for sp in self.spaces]
