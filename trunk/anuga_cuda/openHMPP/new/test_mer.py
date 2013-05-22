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
from hmpp_domain import HMPP_domain

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
mesh_filename = "../merimbula_43200.tsh"   ; x0 = 756000.0 ; x1 = 756500.0
#mesh_filename = "test-100.tsh" ; x0 = 0.25 ; x1 = 0.5
#mesh_filename = "test-20.tsh" ; x0 = 250.0 ; x1 = 350.0
yieldstep = 50
finaltime = 50
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
    

def generate_domain():
    #--------------------------------------------------------------------------
    # Setup Domain only on processor 0
    #--------------------------------------------------------------------------
    domain = create_domain_from_file(mesh_filename, HMPP_domain)
    
    
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


    return domain



#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------

#domain.evolve(yieldstep = yieldstep, finaltime = finaltime)
#t0 = time.time()

#for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
#	domain.write_time()


#barrier()



#--------------------------------------------------
# Merge the individual sww files into one file
#--------------------------------------------------
#domain.sww_merge(delete_old=False)



domain = generate_domain()

from anuga.config import epsilon

if domain.timestepping_method == 'euler':
    timestepping_method = 1
elif domain.timestepping_method == 'rk2':
    timestepping_method = 2
elif domain.timestepping_method == 'rk3':
    timestepping_method = 3
else:
    timestepping_method = 4
print " The timestepping_method is '%s' %d" % (domain.timestepping_method, timestepping_method)


if domain.flow_algorithm == 'tsunami':
    flow_algorithm = 1
elif domain.flow_algorithm == 'yusuke':
    flow_algorithm = 2
else:
    flow_algorithm = 3
print " The flow_algorithm us '%s' %d" % (domain.flow_algorithm, flow_algorithm)


if domain.compute_fluxes_method == 'original':
    compute_fluxes_method = 0
elif domain.compute_fluxes_method == 'wb_1':
    compute_fluxes_method = 1
elif domain.compute_fluxes_method == 'wb_2':
    compute_fluxes_method = 2
elif domain.compute_fluxes_method == 'wb_3':
    compute_fluxes_method = 3
elif domain.compute_fluxes_method == 'tsunami':
    compute_fluxes_method = 4
else:
    compute_fluxes_method = 5
print " The compute_fluxes_method is '%s' %d" % (domain.compute_fluxes_method, compute_fluxes_method)


#import copy
from sort_domain import *

sort_domain(domain)
test_domain = generate_domain()
sort_domain(test_domain)

from hmpp_python_glue import *
from shallow_water_ext import compute_fluxes_ext_central_structure
from shallow_water_ext import gravity_wb as gravity_wb_c
from quantity_ext import extrapolate_second_order_and_limit_by_vertex


domain.convert_boundary_elements()
test_domain.convert_boundary_elements()
                
            

global_cnt = 0
global_correct = True

name_list = ['max_speed', 
                'stage_edge', 'elevation_edge', 'xmom_edge', 'ymom_edge', 
                'height_edge', 'xvelocity_edge', 'yvelocity_edge', 
                
                'stage_centroid', 'elevation_centroid', 'xmom_centroid', 'ymom_centroid', 
                'height_centroid', 'xvelocity_centroid', 'yvelocity_centroid', 
                'friction_centroid',
                
                'stage_vertex', 'elevation_vertex', 'xmom_vertex', 'ymom_vertex', 
                'height_vertex', 'xvelocity_vertex', 'yvelocity_vertex', 
                
                'stage_explicit_update', 'xmom_explicit_update', 'ymom_explicit_update', 
                'stage_semi_implicit_update', 'xmom_semi_implicit_update', 
                'ymom_semi_implicit_update', 

                'stage_x_gradient', 'elevation_x_gradient', 'xmom_x_gradient', 
                'ymom_x_gradient', 'height_x_gradient', 'xvelocity_x_gradient', 
                'yvelocity_x_gradient', 
                
                'stage_y_gradient', 'elevation_y_gradient', 'xmom_y_gradient', 
                'ymom_y_gradient', 'height_y_gradient', 'xvelocity_y_gradient', 
                'yvelocity_y_gradient', 

                'stage_boundary', 'elevation_boundary', 'xmom_boundary', 'ymom_boundary', 
                'height_boundary', 'xvelocity_boundary', 'yvelocity_boundary' 
            ]

def check_result(a, b):
    global global_cnt
    global name_list
    global global_correct
    if not num.allclose(a, b):
        print global_cnt, name_list[global_cnt]
        cnt = 0
        if a.shape.__len__()  == 1:
            for i in range(a.shape[0]):
                if not num.allclose( a[i], b[i] ):
                    cnt += 1
                    if cnt <= 10:
                        print i, a[i], b[i] 
        else:
            for i in range(a.shape[0]):
                if not num.allclose( a[i], b[i] ):
                #if ( a[i] != b[i]).any():
                    cnt += 1
                    if cnt <= 10:
                        print i, a[i], b[i] 
        print cnt
        
        global_correct = False
    global_cnt += 1
            

def check_all(domain, test_domain):
    global global_cnt
    global_cnt = 0
    print domain.flux_timestep, test_domain.flux_timestep
    
    
    check_result(domain.max_speed, test_domain.max_speed)
    
    
    check_result(domain.quantities['stage'].edge_values, 
            test_domain.quantities['stage'].edge_values)
    check_result(domain.quantities['elevation'].edge_values, 
            test_domain.quantities['elevation'].edge_values)
    check_result(domain.quantities['xmomentum'].edge_values, 
            test_domain.quantities['xmomentum'].edge_values)
    check_result(domain.quantities['ymomentum'].edge_values, 
            test_domain.quantities['ymomentum'].edge_values)
    check_result(domain.quantities['height'].edge_values, 
            test_domain.quantities['height'].edge_values)
    check_result(domain.quantities['xvelocity'].edge_values, 
            test_domain.quantities['xvelocity'].edge_values)
    check_result(domain.quantities['yvelocity'].edge_values, 
            test_domain.quantities['yvelocity'].edge_values)
    
    
    check_result(domain.quantities['stage'].centroid_values, 
            test_domain.quantities['stage'].centroid_values)
    check_result(domain.quantities['elevation'].centroid_values, 
            test_domain.quantities['elevation'].centroid_values)
    check_result(domain.quantities['xmomentum'].centroid_values, 
            test_domain.quantities['xmomentum'].centroid_values)
    check_result(domain.quantities['ymomentum'].centroid_values, 
            test_domain.quantities['ymomentum'].centroid_values)
    check_result(domain.quantities['height'].centroid_values, 
            test_domain.quantities['height'].centroid_values)
    check_result(domain.quantities['xvelocity'].centroid_values, 
            test_domain.quantities['xvelocity'].centroid_values)
    check_result(domain.quantities['yvelocity'].centroid_values, 
            test_domain.quantities['yvelocity'].centroid_values)
    check_result(domain.quantities['friction'].centroid_values, 
            test_domain.quantities['friction'].centroid_values)
    
    
    check_result(domain.quantities['stage'].vertex_values, 
            test_domain.quantities['stage'].vertex_values)
    check_result(domain.quantities['elevation'].vertex_values, 
            test_domain.quantities['elevation'].vertex_values)
    check_result(domain.quantities['xmomentum'].vertex_values, 
            test_domain.quantities['xmomentum'].vertex_values)
    check_result(domain.quantities['ymomentum'].vertex_values, 
            test_domain.quantities['ymomentum'].vertex_values)
    check_result(domain.quantities['height'].vertex_values, 
            test_domain.quantities['height'].vertex_values)
    check_result(domain.quantities['xvelocity'].vertex_values, 
            test_domain.quantities['xvelocity'].vertex_values)
    check_result(domain.quantities['yvelocity'].vertex_values, 
            test_domain.quantities['yvelocity'].vertex_values)
    
    
    check_result(domain.quantities['stage'].explicit_update, 
            test_domain.quantities['stage'].explicit_update)
    check_result(domain.quantities['xmomentum'].explicit_update, 
            test_domain.quantities['xmomentum'].explicit_update)
    check_result(domain.quantities['ymomentum'].explicit_update, 
            test_domain.quantities['ymomentum'].explicit_update)
    
        
    check_result(domain.quantities['stage'].semi_implicit_update, 
            test_domain.quantities['stage'].semi_implicit_update)
    check_result(domain.quantities['xmomentum'].semi_implicit_update, 
            test_domain.quantities['xmomentum'].semi_implicit_update)
    check_result(domain.quantities['ymomentum'].semi_implicit_update, 
            test_domain.quantities['ymomentum'].semi_implicit_update)
    
    
    check_result(domain.quantities['stage'].x_gradient, 
            test_domain.quantities['stage'].x_gradient)
    check_result(domain.quantities['elevation'].x_gradient, 
            test_domain.quantities['elevation'].x_gradient)
    check_result(domain.quantities['xmomentum'].x_gradient, 
            test_domain.quantities['xmomentum'].x_gradient)
    check_result(domain.quantities['ymomentum'].x_gradient, 
            test_domain.quantities['ymomentum'].x_gradient)
    check_result(domain.quantities['height'].x_gradient, 
            test_domain.quantities['height'].x_gradient)
    check_result(domain.quantities['xvelocity'].x_gradient, 
            test_domain.quantities['xvelocity'].x_gradient)
    check_result(domain.quantities['yvelocity'].x_gradient, 
            test_domain.quantities['yvelocity'].x_gradient)
    
    
    check_result(domain.quantities['stage'].y_gradient, 
            test_domain.quantities['stage'].y_gradient)
    check_result(domain.quantities['elevation'].y_gradient, 
            test_domain.quantities['elevation'].y_gradient)
    check_result(domain.quantities['xmomentum'].y_gradient, 
            test_domain.quantities['xmomentum'].y_gradient)
    check_result(domain.quantities['ymomentum'].y_gradient, 
            test_domain.quantities['ymomentum'].y_gradient)
    check_result(domain.quantities['height'].y_gradient, 
            test_domain.quantities['height'].y_gradient)
    check_result(domain.quantities['xvelocity'].y_gradient, 
            test_domain.quantities['xvelocity'].y_gradient)
    check_result(domain.quantities['yvelocity'].y_gradient, 
            test_domain.quantities['yvelocity'].y_gradient)
    
    
    check_result(domain.quantities['stage'].boundary_values, 
            test_domain.quantities['stage'].boundary_values)
    check_result(domain.quantities['elevation'].boundary_values, 
            test_domain.quantities['elevation'].boundary_values)
    check_result(domain.quantities['xmomentum'].boundary_values, 
            test_domain.quantities['xmomentum'].boundary_values)
    check_result(domain.quantities['ymomentum'].boundary_values, 
            test_domain.quantities['ymomentum'].boundary_values)
    check_result(domain.quantities['height'].boundary_values, 
            test_domain.quantities['height'].boundary_values)
    check_result(domain.quantities['xvelocity'].boundary_values, 
            test_domain.quantities['xvelocity'].boundary_values)
    check_result(domain.quantities['yvelocity'].boundary_values, 
            test_domain.quantities['yvelocity'].boundary_values)

N = domain.number_of_elements

#for i in range(N):
#    domain.neighbours[i][0] = 1#i+1
#    domain.neighbours[i][1] = 1#i+2
#    domain.neighbours[i][2] = 1#i+3
#
#domain.neighbours[N-3][0] = N-1
#domain.neighbours[N-3][1] = N-2
#domain.neighbours[N-3][2] = N-3
#
#domain.neighbours[N-2][0] = N-1
#domain.neighbours[N-2][1] = N-2
#domain.neighbours[N-2][2] = N-3
#
#domain.neighbours[N-1][0] = N-1
#domain.neighbours[N-1][1] = N-2
#domain.neighbours[N-1][2] = N-3


for i in range(10):
    print "  Initial checking "
    check_all(domain, test_domain)


    #hmpp_distribute_to_vertices_and_edges(
    #hmpp_compute_fluxes(
    hmpp_extrapolate_second_order_and_limit_by_vertex(
            domain,
            yieldstep,
            finaltime,
            0.0,
            epsilon,
            False,
    
            num.int32( compute_fluxes_method ),
            num.int32( flow_algorithm ),
            num.int32( timestepping_method ),
            num.int64( 2),
            num.int32( i)
            )
    
    
    #domain.protect_against_infinitesimal_and_negative_heights()
    #test_domain.distribute_to_vertices_and_edges()
    #test_domain.update_ghosts()
    #test_domain.distribute_to_vertices_and_edges()
    #test_domain.update_boundary()
    #test_domain.update_extrema()
    #test_domain.extrapolate_second_order_sw()
    
    hmpp_extrapolate_second_order_and_limit_by_vertex_normal(
            test_domain,
            yieldstep,
            finaltime,
            0.0,
            epsilon,
            False,
    
            num.int32( compute_fluxes_method ),
            num.int32( flow_algorithm ),
            num.int32( timestepping_method ),
            num.int64( 2),
            num.int32( i)
            )
    #for name in test_domain.conserved_quantities:
    #    Q = test_domain.quantities[name]
    #    #Q.extrapolate_second_order_and_limit_by_vertex()
    #    extrapolate_second_order_and_limit_by_vertex( Q)
    #test_domain.compute_fluxes()
    #compute_fluxes_ext_central_structure(test_domain)
    #gravity_wb_c(test_domain)
    
    
    
    print "\n\n  Final checking  %d \n\n" % i
    check_all(domain, test_domain)

    if not global_correct:
        break
