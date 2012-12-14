

from anuga import Domain
import numpy as num


import pycuda.driver as drv
from pycuda.compiler import SourceModule
auto_init_context = False
if auto_init_context:
    import pycuda.autoinit
else:
    drv.init()
    dev = drv.Device(0)
    ctx = dev.make_context( drv.ctx_flags.MAP_HOST)

from anuga_cuda.config import compute_fluxes_dir

class GPU_domain(Domain):

    def __init__(self, coordinates, vertices,
            boundary=None,
            full_send_dict=None,
            ghost_recv_dict=None,
            number_of_full_nodes=None,
            number_of_full_triangles=None,
            geo_reference=None): 

        Domain.__init__(self,
                coordinates,
                vertices,
                boundary,
                full_send_dict=full_send_dict,
                ghost_recv_dict=ghost_recv_dict,
                number_of_full_nodes=number_of_full_nodes,
                number_of_full_triangles=number_of_full_triangles,
                geo_reference=geo_reference) 

        self.compute_fluxe_mod = SourceModule(
                open("compute_fluxes.cu","r").read(),
                include_dirs=[compute_fluxes_dir]
                )

        self.compute_fluxe_func =mod.get_function(
                "compute_fluxes_central_structure_cuda_single")

        self.gravity_wb_mod = 
        self.gravity_wb_func = mod.get_function("gravity_wb")

    def page_locked_array(self):
        """Convert all necessary array to page-locked array
        """

        self.neighbours = get_page_locked_array(self.neighbours)
        self.neighbour_edges = get_page_locked_array(self.neighbour_edges)
        self.normals = get_page_locked_array(self.normals)
        self.edgelengths = get_page_locked_array(self.edgelengths)
        self.radii = get_page_locked_array(self.radii)
        self.areas = get_page_locked_array(self.areas)
        self.tri_full_flag = get_page_locked_array(self.tri_full_flag)
        self.max_speed = get_page_locked_array(self.max_speed)

        # get edge values
        self.quantities['stage'].edge_values = 
            get_page_locked_array(self.quantities['stage'].edge_values)
        self.quantities['xmomentum'].edge_values = 
            get_page_locked_array(self.quantities['xmomentum'].edge_values)
        self.quantities['ymomentum'].edge_values = 
            get_page_locked_array(self.quantities['ymomentum'].edge_values)
        self.quantities['elevation'].edge_values = 
            get_page_locked_array(self.quantities['elevation'].edge_values)
        
        # get boundary values
        self.quantities['stage'].boundary_values = 
            get_page_locked_array(self.quantities['stage'].boundary_values)
        self.quantities['xmomentum'].boundary_values = 
            get_page_locked_array(
                    self.quantities['xmomentum'].boundary_values)
        self.quantities['ymomentum'].boundary_values = 
            get_page_locked_array(
                    self.quantities['ymomentum'].boundary_values)
        
        # get explicit update
        self.quantities['stage'].explicit_update = 
            get_page_locked_array(self.quantities['stage'].explicit_update)
        self.quantities['xmomentum'].explicit_update = 
            get_page_locked_array(
                    self.quantities['xmomentum'].explicit_update)
        self.quantities['ymomentum'].explicit_update = 
            get_page_locked_array(
                    self.quantities['ymomentum'].explicit_update)
        
        # get vertex values
        self.quantities['stage'].vertex_values = 
            get_page_locked_array(self.quantities['stage'].vertex_values)
        
        # get centroid values
        self.quantities['stage'].centroid_values = 
            get_page_locked_array(self.quantities['stage'].centroid_values)
        self.quantities['elevation'].centroid_values = 
            get_page_locked_array(
                    self.quantities['elevation'].centroid_values)

            
    def evolve(self, 
                yieldstep=None,
                finaltime=None,
                duration=None,
                skip_initial_step=False):

        
        for t in Domain.evolve(self, yieldstep=yieldstep,
                    finaltime=finaltime,
                    duration=duration,
                    skip_initial_step=skip_initial_step):

            yield(t)





def get_page_locked_array(a):
    """This function convert the pageable array
        to page-locked array
    """

    temp_page_lock_p = drv.pagelocked_zeros_like(a,
            mem_flags=drv.host_alloc_flags.DEVICEMAP)
    if len(a.shape) == 1:
        temp_page_lock_p[:] = a
    else:
        temp_page_lock_p[:, :] = a
            
    return temp_page_lock_p
