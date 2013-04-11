# -*- coding: utf-8 -*-

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
                finaltime=0.0,
                duration=0.0,
                skip_initial_step=False):
        
        print self.number_of_elements

        from hmpp_python_glue import hmpp_evolve

        hmpp_evolve(self,
                yieldstep,
                finaltime,
                duration,
                skip_initial_step)



