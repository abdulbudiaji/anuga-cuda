#!/usr/bin/env python

"""Run parallel shallow water domain.

   run using command like:

   mpirun -np m python run_parallel_sw_merimbula.py

   where m is the number of processors to be used.
   
   Will produce sww files with names domain_Pn_m.sww where m is number of processors and
   n in [0, m-1] refers to specific processor that owned this part of the partitioned mesh.
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import os
import sys
import time
import numpy as num

#------------------------
# ANUGA Modules
#------------------------
	
from anuga import Domain
from anuga_cuda import GPU_domain, merimbula_dir
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary

from anuga import rectangular_cross
from anuga import create_domain_from_file




#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

#mesh_filename = "merimbula_10785_1.tsh" ; x0 = 756000.0 ; x1 = 756500.0
mesh_filename = "merimbula_43200.tsh"   ; x0 = 756000.0 ; x1 = 756500.0
#mesh_filename = "test-100.tsh" ; x0 = 0.25 ; x1 = 0.5
#mesh_filename = "test-20.tsh" ; x0 = 250.0 ; x1 = 350.0
mesh_filename = merimbula_dir + mesh_filename
yieldstep = 50
finaltime = 500
verbose = True

#--------------------------------------------------------------------------
# Setup procedures
#--------------------------------------------------------------------------
class Set_Stage:
    """Set an initial condition with constant water height, for x0<x<x1
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))+1.0


class Set_Elevation:
    """Set an elevation
    """

    def __init__(self, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return x/self.h
    

#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
domain = create_domain_from_file(mesh_filename, GPU_domain)

for i in range(len(sys.argv)):
    if sys.argv[i] == '-gpu':
        domain.using_gpu = True
        print " --> Enable GPU version"
    elif sys.argv[i] == '-fs':
        finaltime = float(sys.argv[i+1])
        print " --> Finaltime is reset as %f" % finaltime
    elif sys.argv[i] == '-test':
        domain.cotesting = True
        print " --> Enable Cotesting"
    elif sys.argv[i] == '-ustore':
        domain.store = True
        print " --> Disable storing"

domain.using_gpu = True

domain.set_quantity('stage', Set_Stage(x0, x1, 2.0))
domain.set_datadir('Data')
domain.set_name('merimbula_new')
domain.set_store(True)
    #domain.set_quantity('elevation', Set_Elevation(500.0))

    #print domain.statistics()
    #print domain.get_extent()
    #print domain.get_extent(absolute=True)
    #print domain.geo_reference

#--------------------------------------------------------------------------
# Distribute sequential domain on processor 0 to other processors
#--------------------------------------------------------------------------

#if myid == 0 and verbose: print 'DISTRIBUTING DOMAIN'
#domain = distribute(domain)

#--------------------------------------------------------------------------
# On all processors, setup evolve parameters for domains on all processors
# (all called "domain"
#--------------------------------------------------------------------------

#domain.set_flow_algorithm('2_0')

#domain.smooth = False
#domain.set_default_order(2)
#domain.set_timestepping_method('rk2')
#domain.set_CFL(0.7)
#domain.set_beta(1.5)


#for p in range(numprocs):
#    if myid == p:
#        print 'P%d'%p
#        print domain.get_extent()
#        print domain.get_extent(absolute=True)
#        print domain.geo_reference
#        print domain.s2p_map
#        print domain.p2s_map
#        print domain.tri_l2g
#        print domain.node_l2g
#    else:
#        pass
#
#    barrier()


#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
#------------------------------------------------------------------------------
Br = Reflective_boundary(domain)      # Solid reflective wall

domain.set_boundary({'outflow' :Br, 'inflow' :Br, 'inner' :Br, 'exterior' :Br, 'open' :Br})

#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------

t0 = time.time()

for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
	domain.write_time()


#barrier()



#--------------------------------------------------
# Merge the individual sww files into one file
#--------------------------------------------------
#domain.sww_merge(delete_old=False)





