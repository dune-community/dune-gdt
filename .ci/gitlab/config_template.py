#!/usr/bin/env python3
# kate: indent-width 2;

import os
import jinja2
from itertools import product
import urllib.request

template_url = 'https://raw.githubusercontent.com/dune-community/dune-xt-ci/master/xt_and_gdt_config_template.txt'
tpl = urllib.request.urlopen(template_url).read().decode('utf-8').replace('DUNE_XT_OR_DUNE_GDT', 'dune-gdt')
tpl = jinja2.Template(tpl)
images = ['debian-unstable_gcc_full', 'debian_gcc_full', 'debian_clang_full']
subdirs = ['gdt']
kinds = ['cpp', 'headercheck']
matrix = product(images, subdirs, kinds)
with open(os.path.join(os.path.dirname(__file__), 'config.yml'), 'wt') as yml:
    yml.write(tpl.render(matrix=matrix, images=images, kinds=kinds, subdirs=subdirs))
