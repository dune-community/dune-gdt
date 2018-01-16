from __future__ import division

tpl_str='''#!/bin/bash
##
## optional: energy policy tags
#@ energy_policy_tag = NONE
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = {{ 'test' if nodes < 17 else 'general' }}
#@ node = {{ nodes }}
#@ island_count=1,1
# other example
#@ tasks_per_node = 28
#@ wall_clock_limit = 0:12:30
##                    1 h 20 min 30 secs
#@ job_name = speedup_{{ MACRO }}_{{ nodes }}N_1C
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/build/main_gdt/gcc/dune-gdt/dune/gdt/test/
#@ output = speedup_blockswipdg_{{ MACRO }}x__{{ nodes }}N_1C_$(jobid).out
#@ error = speedup_blockswipdg_{{ MACRO }}x__{{ nodes }}N_1C_$(jobid).err
#@ notification=always
#@ notify_user=rene.milk@wwu.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
source $HOME/.modules
#optional: 
#module load mpi_pinning/hybrid_blocked

BIN=$HOME/build/main_gdt/gcc/dune-gdt/dune/gdt/test/test_block_swipdg_discretization


OPT="/home/hpc/pr62zo/di73dez2/main_gdt/dune-gdt/dune/gdt/test/block_swipdg_discretization.ini \
-global.datadir $HOME/results/block-yaspgrid_speedup_n{{ procs }}_{{ MACRO }}x__T1 "

mpirun -n {{ procs }} $BIN ${OPT}
ERR=$(ll speedup_blockswipdg_{{ MACRO }}x__{{ nodes }}N_1C_$(jobid).err)
push_notify "${ERR} {{ procs }}"
'''

from jinja2 import Template
tpl=Template(tpl_str)

MACRO=64
def output_run(i):
  nodes = int(pow(2,i))
  procs = 28 * nodes
  #print('nodes {} - procs {}'.format(nodes, procs))
  fn = 'node_{}.submit'.format(nodes)
  open(fn, 'wt').write(tpl.render(**locals()))
  print('llsubmit {}'.format(fn))


for i in range(1,3):
  output_run(i)

