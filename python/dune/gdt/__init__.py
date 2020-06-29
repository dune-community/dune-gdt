# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017 - 2018)
#   René Fritze     (2016, 2018)
# ~~~

from dune.xt import guarded_import

for mod_name in ( # order should not matter!
        '_discretefunction_discretefunction',
        '_discretefunction_dof_vector',
        '_functionals_interfaces_common',
        '_functionals_interfaces_eigen',
        '_functionals_interfaces_istl',
        '_functionals_vector_based',
        '_local_bilinear_forms_element_integrals',
        '_local_bilinear_forms_element_interface',
        '_local_functionals_element_integrals',
        '_local_functionals_element_interface',
        '_local_functionals_intersection_integrals',
        '_local_functionals_intersection_interface',
        '_local_integrands_advection',
        '_local_integrands_binary_element_interface',
        '_local_integrands_element_product',
        '_local_integrands_intersection_product',
        '_local_integrands_laplace',
        '_local_integrands_unary_element_interface',
        '_operators_interfaces_common',
        '_operators_interfaces_eigen',
        '_operators_interfaces_istl',
        '_operators_matrix_based',
        '_spaces_h1_continuous_lagrange',
        '_spaces_interface',
        '_spaces_l2_discontinuous_lagrange',
        '_tools_dirichlet_constraints',
        '_tools_sparsity_pattern',
):
    guarded_import(globals(), 'dune.gdt', mod_name)
