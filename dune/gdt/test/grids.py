class Grids(object):
    def __init__(self, cache):
        try:
            have_alugrid = cache['dune-alugrid']
        except KeyError:
            have_alugrid = False
        self.all_parts_fmt = [self.alu_conf_part_fmt, self.alu_cube_part_fmt, self.yasp_part_fmt] if have_alugrid else [
            self.yasp_part_fmt]
        self.all_views_fmt = [self.alu_conf_view_fmt, self.alu_cube_view_fmt, self.yasp_view_fmt] if have_alugrid else [
            self.yasp_view_fmt]


class LeafGrids(Grids):
    world_dim = [2, 3]

    alu_cube_part_fmt = 'AluCube{}dLeafGridPartType'
    alu_conf_part_fmt = 'AluConform{}dLeafGridPartType'
    yasp_part_fmt = 'Yasp{}dLeafGridPartType'

    alu_cube_view_fmt = 'AluCube{}dLeafGridViewType'
    alu_conf_view_fmt = 'AluConform{}dLeafGridViewType'
    yasp_view_fmt = 'Yasp{}dLeafGridViewType'

    def __init__(self, cache):
        super(LeafGrids, self).__init__(cache)


class LevelGrids(Grids):
    world_dim = [2, 3]

    alu_cube_part_fmt = 'AluCube{}dLevelGridPartType'
    alu_conf_part_fmt = 'AluConform{}dLevelGridPartType'
    yasp_part_fmt = 'Yasp{}dLevelGridPartType'

    alu_cube_view_fmt = 'AluCube{}dLevelGridViewType'
    alu_conf_view_fmt = 'AluConform{}dLevelGridViewType'
    yasp_view_fmt = 'Yasp{}dLevelGridViewType'

    def __init__(self, cache):
        super(LevelGrids, self).__init__(cache)