# -*- coding: utf-8 -*-


import ctypes
import sys


import numpy


# This is used to fix the dlopen() issues for python
# and make the reference links (like OpenHMPP global)
flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

from hmpp_python_glue import hmpp_python_test
sys.setdlopenflags(flags)


from anuga import Domain

import numpy as num



class HMPP_domain(Domain):

    def __init__(self, 
                coordinates=None, 
                vertices=None,
                 boundary=None,
                 source=None,
                 triangular=None,
                 conserved_quantities=None,
                 evolved_quantities=None,
                 other_quantities=None,
                 tagged_elements=None,
                 geo_reference=None,
                 use_inscribed_circle=False,
                 mesh_fulename=None,
                 use_cache=False,
                 verbose=False,
                 full_send_dict=None,
                 ghost_recv_dict=None,
                 starttime=0.0,
                 processor=0,
                 numproc=1,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None,
                 ghost_layer_width=2
                 ): #jj added this

        Domain.__init__(self,
                        coordinates,
                        vertices,
                        boundary,
                        full_send_dict=full_send_dict,
                        ghost_recv_dict=ghost_recv_dict,
                        number_of_full_nodes=number_of_full_nodes,
                        number_of_full_triangles=number_of_full_triangles,
                        geo_reference=geo_reference) #jj added this



    def evolve(self,
                yieldstep=0.0,
                finaltime=1.0,
                duration=0.0,
                skip_initial_step=False):
        

        if self.store is True and self.get_time() == self.get_starttime():
            self.initialise_storage()

        from hmpp_python_glue import hmpp_evolve
        from anuga.config import epsilon

        if self.timestepping_method == 'euler':
            timestepping_method = 1
        elif self.timestepping_method == 'rk2':
            timestepping_method = 2
        elif self.timestepping_method == 'rk3':
            timestepping_method = 3
        else:
            timestepping_method = 4
        print " The timestepping_method is '%s' %d" % (self.timestepping_method, timestepping_method)


        if self.flow_algorithm == 'tsunami':
            flow_algorithm = 1
        elif self.flow_algorithm == 'yusuke':
            flow_algorithm = 2
        else:
            flow_algorithm = 3
        print " The flow_algorithm us '%s' %d" % (self.flow_algorithm, flow_algorithm)
       

        if self.compute_fluxes_method == 'original':
            compute_fluxes_method = 0
        elif self.compute_fluxes_method == 'wb_1':
            compute_fluxes_method = 1
        elif self.compute_fluxes_method == 'wb_2':
            compute_fluxes_method = 2
        elif self.compute_fluxes_method == 'wb_3':
            compute_fluxes_method = 3
        elif self.compute_fluxes_method == 'tsunami':
            compute_fluxes_method = 4
        else:
            compute_fluxes_method = 5
        print " The compute_fluxes_method is '%s' %d" % (self.compute_fluxes_method, compute_fluxes_method)

            
        yield_step = 0
        while True :
            tmp_timestep = hmpp_evolve(self,
                    yieldstep,
                    finaltime,
                    duration,
                    epsilon,
                    skip_initial_step,

                    numpy.int32( compute_fluxes_method ),
                    numpy.int32( flow_algorithm ),
                    numpy.int32( timestepping_method ),
                    numpy.int32( yield_step)
                    )

            yield_step = 1
            print " Python: tmp_timestep %lf " % tmp_timestep

            cmd = raw_input("HMPP_DOMAIN: Quit?[q]")
            if 'q' in cmd or 'Q' in cmd:
                break
            if tmp_timestep >= finaltime - epsilon: 
                print " Evolve finish"
                break

