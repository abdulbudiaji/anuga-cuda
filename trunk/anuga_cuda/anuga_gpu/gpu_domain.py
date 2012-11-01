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





    def compute_fluxes(self):
        """Compute fluxes and timestep suitable for all volumes in domain.
 
        Compute total flux for each conserved quantity using "flux_function"
 
        Fluxes across each edge are scaled by edgelengths and summed up
        Resulting flux is then scaled by area and stored in
        explicit_update for each of the three conserved quantities
        stage, xmomentum and ymomentum
 
        The maximal allowable speed computed by the flux_function for each volume
        is converted to a timestep that must not be exceeded. The minimum of
        those is computed as the next overall timestep.
 
        Post conditions:
          domain.explicit_update is reset to computed flux values
          domain.timestep is set to the largest step satisfying all volumes.
 
        This wrapper calls the underlying C version of compute fluxes
        """
 
        import sys
        from gpu_python_glue import compute_fluxes_ext_central_new_gpu as compute_fluxes_ext
 
        # Shortcuts
        Stage = self.quantities['stage']
        Xmom = self.quantities['xmomentum']
        Ymom = self.quantities['ymomentum']
        Bed = self.quantities['elevation']
 
        timestep = float(sys.maxint)
        print timestep, self, Stage, Xmom, Ymom, Bed
        print self.tri_full_flag
        print "areas: ", self.areas

        flux_timestep = compute_fluxes_ext(timestep, self, Stage, Xmom, Ymom, Bed)
        self.flux_timestep = flux_timestep
        print "Updates: "
        print "stage ", Stage.explicit_update
        print "xmom  ", Xmom.explicit_update
        print "ymom  ", Ymom.explicit_update

