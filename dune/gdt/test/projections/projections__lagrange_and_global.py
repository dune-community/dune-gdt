from dune.xt.codegen import typeid_to_typedef_name

import spaces as sp

cg = sp.CG(cache)
spaces = cg.spaces
names = cg.names
spaces_with_names = zip(spaces, names)
