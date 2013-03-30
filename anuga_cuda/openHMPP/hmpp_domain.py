# -*- coding: utf-8 -*-
"""Class Parallel_shallow_water_domain -
2D triangular domains for finite-volume computations of
the shallow water equation, with extra structures to allow
communication between other Parallel_domains and itself

This module contains a specialisation of class Domain
from module shallow_water.py

Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2004-2005

"""

from anuga import Domain

import numpy as num



class GPU_domain(Domain):

    def __init__(self, coordinates, vertices,
                 boundary=None,
                 full_send_dict=None,
                 ghost_recv_dict=None,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None,
                 geo_reference=None): #jj added this

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

        from hmpp_python_glue import hmpp_evolve

        hmpp_evolve(self,
                yieldstep,
                finaltime,
                duration)


