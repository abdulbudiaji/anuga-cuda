#!/usr/bin/env python

from anuga import Domain
import numpy 


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
from anuga_cuda.config import gravity_dir
from anuga_cuda.config import extrapolate_dir
from anuga_cuda.config import protect_dir
from anuga_cuda.config import balance_dir
from anuga_cuda.config import interpolate_dir

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

        from anuga_cuda.merimbula_data.sort_domain import sort_domain
        sort_domain(self)

        self.using_gpu = True

        self.end_event = drv.Event()

        self.timestep_array = drv.pagelocked_zeros(self.number_of_elements,
            dtype = numpy.float64,
            mem_flags=drv.host_alloc_flags.DEVICEMAP)
       
        # compute_fluxes function
        self.compute_fluxes_mod = SourceModule(
                open(compute_fluxes_dir+"compute_fluxes.cu","r").read(),
                include_dirs=[compute_fluxes_dir]
                )

        self.compute_fluxes_func = self.compute_fluxes_mod.get_function(
                "compute_fluxes_central_structure_cuda_single")

        # gravity_wb function
        self.gravity_wb_mod = SourceModule(
                open(gravity_dir+"gravity.cu","r").read(),
                include_dirs=[gravity_dir]
                )

        self.gravity_wb_func = self.gravity_wb_mod.get_function(
                "gravity_wb")

        # extrapolate_second_order_sw function
        self.extrapolate_second_order_sw_mod = SourceModule(
                open(extrapolate_dir+"extrapolate_second_order_sw.cu").read(),
                options=["--ptxas-options=-v"],
                include_dirs=[extrapolate_dir]
            )

        self.extrapolate_second_order_sw_true_func = \
            self.extrapolate_second_order_sw_mod.get_function(
                    "extrapolate_second_order_sw_true")

        self.extrapolate_second_order_sw_false_func = \
            self.extrapolate_second_order_sw_mod.get_function(
                "extrapolate_second_order_sw_false")

        # extrapolate_second_order_and_limit_by_vertex function
        self.extrapolate_second_order_and_limit_by_vertex_mod = \
            SourceModule(
                open(extrapolate_dir + \
                    "extrapolate_second_order_and_limit_by_vertex.cu"
                    ).read(),
                include_dirs=[extrapolate_dir]
                )

        self.extrapolate_second_order_and_limit_by_vertex_func = \
            self.extrapolate_second_order_and_limit_by_vertex_mod.get_function(
                    "extrapolate_second_order_and_limit_by_vertex")

        # protect function
        self.protect_mod = SourceModule(
                open(protect_dir + "protect.cu").read(),
                include_dirs=[protect_dir]
                )

        self.protect_sw_func = self.protect_mod.get_function(
                "_protect_sw")
        self.protect_swb2_func = self.protect_mod.get_function(
                "_protect_swb2")
        
        # balance function
        self.balance_deep_and_shallow_mod = SourceModule(
                open(balance_dir + "balance_deep_and_shallow.cu").read(),
                include_dirs=[balance_dir]
                )

        self.balance_deep_and_shallow_func = \
            self.balance_deep_and_shallow_mod.get_function(
                "_balance_deep_and_shallow")


    def lock_array_page(self):
        self.neighbours = get_page_locked_array(self.neighbours)
        self.neighbour_edges = get_page_locked_array(self.neighbour_edges)
        self.normals = get_page_locked_array(self.normals)
        self.edgelengths = get_page_locked_array(self.edgelengths)
        self.radii = get_page_locked_array(self.radii)
        self.areas = get_page_locked_array(self.areas)
        self.tri_full_flag = get_page_locked_array(self.tri_full_flag)
        self.max_speed = get_page_locked_array(self.max_speed)
        self.vertex_coordinates = get_page_locked_array(self.vertex_coordinates)

        # get edge values
        self.quantities['stage'].edge_values = \
            get_page_locked_array(self.quantities['stage'].edge_values)
        self.quantities['xmomentum'].edge_values = \
            get_page_locked_array(self.quantities['xmomentum'].edge_values)
        self.quantities['ymomentum'].edge_values = \
            get_page_locked_array(self.quantities['ymomentum'].edge_values)
        self.quantities['elevation'].edge_values = \
            get_page_locked_array(self.quantities['elevation'].edge_values)
        
        # get boundary values
        self.quantities['stage'].boundary_values = \
            get_page_locked_array(self.quantities['stage'].boundary_values)
        self.quantities['xmomentum'].boundary_values = \
            get_page_locked_array(
                    self.quantities['xmomentum'].boundary_values)
        self.quantities['ymomentum'].boundary_values = \
            get_page_locked_array(
                    self.quantities['ymomentum'].boundary_values)
        
        # get explicit update
        self.quantities['stage'].explicit_update = \
            get_page_locked_array(self.quantities['stage'].explicit_update)
        self.quantities['xmomentum'].explicit_update = \
            get_page_locked_array(
                    self.quantities['xmomentum'].explicit_update)
        self.quantities['ymomentum'].explicit_update = \
            get_page_locked_array(
                    self.quantities['ymomentum'].explicit_update)
        
        # get vertex values
        self.quantities['stage'].vertex_values = \
            get_page_locked_array(self.quantities['stage'].vertex_values)
        
        # get centroid values
        self.quantities['stage'].centroid_values = \
            get_page_locked_array(self.quantities['stage'].centroid_values)
        self.quantities['elevation'].centroid_values = \
            get_page_locked_array(
                    self.quantities['elevation'].centroid_values)
    
    def allocate_device_array(self):
        """Convert all necessary array to page-locked array
        """

        # auxiliary arrays
        self.timestep_array_gpu = get_device_array(self.timestep_array)
        self.stage_centroid_store_gpu = get_device_array(
                self.quantities['stage'].centroid_values)
        self.xmomentum_centroid_store_gpu = get_device_array(
                self.quantities['xmomentum'].centroid_values)
        self.ymomentum_centroid_store_gpu = get_device_array(
                self.quantities['ymomentum'].centroid_values)


        self.neighbours_gpu = get_device_array(self.neighbours)
        self.neighbour_edges_gpu = get_device_array(self.neighbour_edges)
        self.surrogate_neighbours_gpu = \
            get_device_array(self.surrogate_neighbours)
        self.normals_gpu = get_device_array(self.normals)
        self.edgelengths_gpu = get_device_array(self.edgelengths)
        self.radii_gpu = get_device_array(self.radii)
        self.areas_gpu = get_device_array(self.areas)
        self.tri_full_flag_gpu = get_device_array(self.tri_full_flag)
        self.max_speed_gpu = get_device_array(self.max_speed)
        
        # domain coordinates
        self.vertex_coordinates_gpu = get_device_array(self.vertex_coordinates)
        self.edge_coordinates_gpu = get_device_array(self.edge_coordinates)
        self.centroid_coordinates_gpu = \
            get_device_array(self.centroid_coordinates)
        
        self.number_of_boundaries_gpu = \
            get_device_array(self.number_of_boundaries)

        # get edge values
        self.quantities['stage'].edge_values_gpu = \
            get_device_array(self.quantities['stage'].edge_values)
        self.quantities['xmomentum'].edge_values_gpu = \
            get_device_array(self.quantities['xmomentum'].edge_values)
        self.quantities['ymomentum'].edge_values_gpu = \
            get_device_array(self.quantities['ymomentum'].edge_values)
        self.quantities['elevation'].edge_values_gpu = \
            get_device_array(self.quantities['elevation'].edge_values)
        
        # get boundary values
        self.quantities['stage'].boundary_values_gpu = \
            get_device_array(self.quantities['stage'].boundary_values)
        self.quantities['xmomentum'].boundary_values_gpu = \
            get_device_array(
                    self.quantities['xmomentum'].boundary_values)
        self.quantities['ymomentum'].boundary_values_gpu = \
            get_device_array(
                    self.quantities['ymomentum'].boundary_values)
        
        # get explicit update
        self.quantities['stage'].explicit_update_gpu = \
            get_device_array(self.quantities['stage'].explicit_update)
        self.quantities['xmomentum'].explicit_update_gpu = \
            get_device_array(
                    self.quantities['xmomentum'].explicit_update)
        self.quantities['ymomentum'].explicit_update_gpu = \
            get_device_array(
                    self.quantities['ymomentum'].explicit_update)
        
        # get vertex values
        self.quantities['stage'].vertex_values_gpu = \
            get_device_array(self.quantities['stage'].vertex_values)
        self.quantities['elevation'].vertex_values_gpu = \
            get_device_array(
                    self.quantities['elevation'].vertex_values)
        self.quantities['xmomentum'].vertex_values_gpu = \
            get_device_array(self.quantities['xmomentum'].vertex_values)
        self.quantities['ymomentum'].vertex_values_gpu = \
            get_device_array(self.quantities['ymomentum'].vertex_values)
        
        # get centroid values
        self.quantities['stage'].centroid_values_gpu = \
            get_device_array(self.quantities['stage'].centroid_values)
        self.quantities['elevation'].centroid_values_gpu = \
            get_device_array(
                    self.quantities['elevation'].centroid_values)
        self.quantities['xmomentum'].centroid_values_gpu = \
            get_device_array(self.quantities['xmomentum'].centroid_values)
        self.quantities['ymomentum'].centroid_values_gpu = \
            get_device_array(self.quantities['ymomentum'].centroid_values)
        
        for name in self.conserved_quantities:
            """ ['stage', 'xmomentum', 'ymomentum']"""
            Q = self.quantities[name]
            Q.x_gradient_gpu = get_de(Q.x_gradient)
            Q.y_gradient_gpu = get_de(Q.y_gradient)


    def asynchronous_transfer(self):
        # auxiliary arrays
        asy_cpy(self.timestep_array, self.timestep_array_gpu)

        asy_cpy(self.neighbours, self.neighbours_gpu)
        asy_cpy(self.neighbour_edges, self.neighbour_edges_gpu)
        asy_cpy(self.surrogate_neighbours, self.surrogate_neighbours_gpu)
        asy_cpy(self.normals, self.normals_gpu)
        asy_cpy(self.edgelengths, self.edgelengths_gpu)
        asy_cpy(self.radii, self.radii_gpu)
        asy_cpy(self.areas, self.areas_gpu)
        asy_cpy(self.tri_full_flag, self.tri_full_flag_gpu)
        asy_cpy(self.max_speed, self.max_speed_gpu)

        # domain coordinates
        asy_cpy(self.vertex_coordinates, self.vertex_coordinates_gpu)
        asy_cpy(self.edge_coordinates, self.edge_coordinates_gpu)
        asy_cpy(self.centroid_coordinates, self.centroid_coordinates_gpu)

        asy_cpy(self.number_of_boundaries, self.number_of_boundaries_gpu)

        # get edge values
        asy_cpy(
                self.quantities['stage'].edge_values,
                self.quantities['stage'].edge_values_gpu)
        asy_cpy(
                self.quantities['xmomentum'].edge_values,
                self.quantities['xmomentum'].edge_values_gpu)
        asy_cpy(
                self.quantities['ymomentum'].edge_values,
                self.quantities['ymomentum'].edge_values_gpu)
        asy_cpy(
                self.quantities['elevation'].edge_values,
                self.quantities['elevation'].edge_values_gpu)

        # get boundary values
        asy_cpy(
                self.quantities['stage'].boundary_values,
                self.quantities['stage'].boundary_values_gpu)
        asy_cpy(
                self.quantities['xmomentum'].boundary_values,
                self.quantities['xmomentum'].boundary_values_gpu)
        asy_cpy(
                self.quantities['ymomentum'].boundary_values,
                self.quantities['ymomentum'].boundary_values_gpu)

        # get explicit update
        asy_cpy(
                self.quantities['stage'].explicit_update,
                self.quantities['stage'].explicit_update_gpu)
        asy_cpy(
                self.quantities['xmomentum'].explicit_update,
                self.quantities['xmomentum'].explicit_update_gpu)
        asy_cpy(
                self.quantities['ymomentum'].explicit_update,
                self.quantities['ymomentum'].explicit_update_gpu)

        # get vertex values
        asy_cpy(
                self.quantities['stage'].vertex_values,
                self.quantities['stage'].vertex_values_gpu)
        asy_cpy(
                self.quantities['elevation'].vertex_values,
                self.quantities['elevation'].vertex_values_gpu)
        asy_cpy(
                self.quantities['xmomentum'].vertex_values,
                self.quantities['xmomentum'].vertex_values_gpu)
        asy_cpy(
                self.quantities['ymomentum'].vertex_values,
                self.quantities['ymomentum'].vertex_values_gpu)

        # get centroid values
        asy_cpy(
                self.quantities['stage'].centroid_values,
                self.quantities['stage'].centroid_values_gpu)
        asy_cpy(
                self.quantities['elevation'].centroid_values,
                self.quantities['elevation'].centroid_values_gpu)
        asy_cpy(
                self.quantities['xmomentum'].centroid_values,
                self.quantities['xmomentum'].centroid_values_gpu)
        asy_cpy(
                self.quantities['ymomentum'].centroid_values,
                self.quantities['ymomentum'].centroid_values_gpu)





    def compute_fluxes(self):
        if self.using_gpu :
            W1 = 32
            W2 = 1
            W3 = 1

            strm_cf = drv.Stream()

            self.compute_fluxes_func(
                numpy.uint64(self.number_of_elements),
                numpy.float64(self.g),
                numpy.float64(self.epsilon),
                numpy.float64(self.H0 * self.H0),
                numpy.float64(self.H0 * 10),
                numpy.uint32(self.optimise_dry_cells),
                self.timestep_array_gpu,
                self.neighbours_gpu,
                self.neighbour_edges_gpu,
                self.normals_gpu,
                self.edgelengths_gpu,
                self.radii_gpu,
                self.areas_gpu,
                self.tri_full_flag_gpu,
                self.quantities['stage'].edge_values_gpu,
                self.quantities['xmomentum'].edge_values_gpu,
                self.quantities['ymomentum'].edge_values_gpu,
                self.quantities['elevation'].edge_values_gpu,
                self.quantities['stage'].boundary_values_gpu,
                self.quantities['xmomentum'].boundary_values_gpu,
                self.quantities['ymomentum'].boundary_values_gpu,
                self.quantities['stage'].explicit_update_gpu,
                self.quantities['xmomentum'].explicit_update_gpu,
                self.quantities['ymomentum'].explicit_update_gpu,
                self.max_speed_gpu,
                block = (W1, W2, W3),
                grid = ((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3), 1),
                stream = strm_cf)
                
            strm_g = drv.Stream()

            self.gravity_wb_func(
                self.quantities['stage'].vertex_values_gpu,
                self.quantities['stage'].edge_values_gpu,
                self.quantities['stage'].centroid_values_gpu,
                self.quantities['elevation'].edge_values_gpu,
                self.quantities['elevation'].centroid_values_gpu,
                self.vertex_coordinates_gpu,
                self.quantities['xmomentum'].explicit_update_gpu,
                self.quantities['ymomentum'].explicit_update_gpu,
                self.normals_gpu,
                self.areas_gpu,
                self.edgelengths_gpu,
                numpy.float64(self.g),
                numpy.uint64(self.number_of_elements),
                block = (W1, W2, W3),
                grid =((self.number_of_elements + W1*W2*W3 -1)/(W1*W2*W3),1),
                stream =strm_g)

            drv.memcpy_dtoh_async(
                    self.timestep_array, 
                    self.timestep_array_gpu, 
                    strm_cf)

            strm_cf.synchronize()

            b = numpy.argsort( self.timestep_array)
            self.flux_timestep = self.timestep_array[ b[0] ]
        else:
            Domain.compute_fluxes(self)
           
    def balance_deep_and_shallow(self):
        if  self.using_gpu:
            W1 = 32
            W2 = 1
            W3 = 1
            self.balance_deep_and_shallow_func(
                numpy.float64(self.H0),
                numpy.float64(self.alpha_balance),

                numpy.int32(self.tight_slope_limiters),
                numpy.int32(self.use_centroid_velocities),
                
                self.quantities['stage'].centroid_values_gpu,
                self.quantities['elevation'].centroid_values_gpu,
                self.quantities['xmomentum'].centroid_values_gpu,
                self.quantities['ymomentum'].centroid_values_gpu,
                
                self.quantities['stage'].vertex_values_gpu,
                self.quantities['elevation'].vertex_values_gpu,
                self.quantities['xmomentum'].vertex_values_gpu,
                self.quantities['ymomentum'].vertex_values_gpu,

                block = (W1, W2, W3),
                grid = ((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3), 1),
                )
        else:
            Domain.balance_deep_and_shallow(self)

    def distribute_to_vertices_and_edges(self):
        if  self.using_gpu:
            W1 = 32
            W2 = 1
            W3 = 1
            if self.compute_fluxes_method == 'tsunami':
                self.protect_swb2_func(
                        numpy.int32(self.number_of_elements),
                        numpy.float64(self.minimum_allowed_height),
                        numpy.float64(self.maximum_allowed_speed),
                        number_of_elements.float64(self.epsilon),
                        self.quantities['stage'].centroid_values_gpu,
                        self.quantities['stage'].vertex_values_gpu,
                        self.quantities['elevation'].centroid_values_gpu,
                        self.quantities['elevation'].vertex_values_gpu,
                        self.quantities['xmomentum'].centroid_values_gpu,
                        self.quantities['ymomentum'].centroid_values_gpu,
                        self.areas_gpu,
                        block = (W1, W2, W3),
                        grid=((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                        )

                self.extrapolate_second_order_edge_swb2_func(
                        numpy.int32(self.number_of_elements),
                        numpy.int32(self.optimise_dry_cells),
                        numpy.int32(self.extrapolate_velocity_second_order),
                        numpy.float64(self.minimum_allowed_height),
                        numpy.float64(self.beta_w),
                        numpy.float64(self.beta_w_dry),
                        numpy.float64(self.beta_uh),
                        numpy.float64(self.beta_uh_dry),
                        numpy.float64(self.beta_vh),
                        numpy.float64(self.beta_vh_dry),

                        self.surrogate_neighbours_gpu,
                        self.number_of_boundaries_gpu,
                        self.centroid_coordinates_gpu,
                        self.quantities['stage'].centroid_values_gpu,
                        self.quantities['elevation'].centroid_values_gpu,
                        self.quantities['xmomentum'].centroid_values_gpu,
                        self.quantities['ymomentum'].centroid_values_gpu,
                        
                        self.edge_coordinates_gpu,
                        
                        self.quantities['stage'].edge_values_gpu,
                        self.quantities['elevation'].edge_values_gpu,
                        self.quantities['xmomentum'].edge_values_gpu,
                        self.quantities['ymomentum'].edge_values_gpu,
                        
                        self.quantities['stage'].vertex_values_gpu,
                        self.quantities['elevation'].vertex_values_gpu,
                        self.quantities['xmomentum'].vertex_values_gpu,
                        self.quantities['ymomentum'].vertex_values_gpu,
                        self.stage_centroid_store_gpu,
                        self.xmomentum_centroid_store_gpu,
                        self.ymomentum_centroid_store_gpu,
                        block = (W1, W2, W3),
                        grid=((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                        )
                
            elif self.use_edge_limiter:
                self.protect_against_infinitesimal_and_negative_heights()

                for name in self.conserved_quantities:
                    Q = self.quantities[name]
                    if self._order_ == 1:
                        #Q.extrapolate_first_order()
                        self.extrapolate_first_order_func(
                            Q.centroid_values_gpu,
                            Q.edge_values_gpu,
                            Q.vertex_values_gpu,
                            block = (W1, W2, W3),
                            grid =((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                            )

                        drv.memset_d32(Q.x_gradient_gpu,0,self.number_of_elements*2)
                        drv.memset_d32(Q.y_gradient_gpu,0,self.number_of_elements*2)
                        

                    elif self._order_ == 2:
                        #Q.extrapolate_second_order_and_limit_by_edge()
                    else:
                        raise Exception('Unknown order')


                self.balance_deep_and_shallow()

                for name in self.conserved_quantities:
                    Q = self.quantities[name]
                    #Q.interpolate_from_vertices_to_edges()
                    self.interpolate_from_vertices_to_edges_func(
                            numpy.int32(Q.vertex_values.shape[0]),
                            Q.vertex_values_gpu,
                            Q.edge_values_gpu,
                            block = (W1, W2, W3),
                            grid =((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                            )



            # using vertex limiter
            else:
                self.protect_against_infinitesimal_and_negative_heights()

                if self.optimised_gradient_limiter:
                    if self._order_ == 1:
                        for name in self.conserved_quantities:
                            Q = self.quantities[name]
                            Q.extrapolate_first_order_sw()
                    elif self._order_ == 2:
                        self.extrapolate_second_order_sw()
                    else:
                        raise Exception('Unknown order')

                else:
                    for name in self.conserved_quantities:
                        Q = self.quantities[name]

                        if self._order_ == 1:
                            #Q.extrapolate_first_order()
                            self.extrapolate_first_order_func(
                                Q.centroid_values_gpu,
                                Q.edge_values_gpu,
                                Q.vertex_values_gpu,
                                block = (W1, W2, W3),
                                grid =(
                                    (self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                                  )

                            drv.memset_d32(
                                Q.x_gradient_gpu,0,self.number_of_elements*2)
                            drv.memset_d32(
                                Q.y_gradient_gpu,0,self.number_of_elements*2)
                        elif self._order_ == 2:
                            Q.extrapolate_second_order_and_limit_by_vertex()
                            #FIXME
                            #self.extrapolate_second_order_and_limit_by_vertex_func(
                            #    numpy.float64( Q.beta),
                            #    self.centroid_coordinates_gpu,
                            #    self.vertex_coordinates_gpu,
                            #    self.number_of_boundaries_gpu,
                            #    self.surrogate_neighbours_gpu,
                            #    self.neighbours_gpu,


                        else:
                            raise Exception('Unknown order')

                self.balance_deep_and_shallow()

                for name in self.conserved_quantities:
                    Q = self.quantities[name]
                    #Q.interpolate_from_vertices_to_edges()
                    self.interpolate_from_vertices_to_edges_func(
                            numpy.int32(Q.vertex_values.shape[0]),
                            Q.vertex_values_gpu,
                            Q.edge_values_gpu,
                            block = (W1, W2, W3),
                            grid =((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                            )
        else:
            Domain.distribute_to_vertices_and_edges(self)
            

    def protect_against_infinitesimal_and_negative_heights(self):
        if  self.using_gpu:
            W1 = 32
            W2 = 1
            W3 = 1

            if self.flow_algorithm == 'tsunami':
                self.protect_swb2_func(
                        numpy.int32(self.number_of_elements),
                        numpy.float64(self.minimum_allowed_height),
                        numpy.float64(self.maximum_allowed_speed),
                        number_of_elements.float64(self.epsilon),
                        self.quantities['stage'].centroid_values_gpu,
                        self.quantities['stage'].vertex_values_gpu,
                        self.quantities['elevation'].centroid_values_gpu,
                        self.quantities['elevation'].vertex_values_gpu,
                        self.quantities['xmomentum'].centroid_values_gpu,
                        self.quantities['ymomentum'].centroid_values_gpu,
                        self.areas_gpu,
                        block = (W1, W2, W3),
                        grid=((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                        )
                    
            else:
                self.protect_sw_func(
                        numpy.float64( self.minimum_allowed_height), 
                        numpy.float64( self.maximum_allowed_speed),
                        numpy.float64( self.epsilon), 
                        self.quantities['stage'].centroid_values_gpu,
                        self.quantities['elevation'].centroid_values_gpu, 
                        self.quantities['xmomentum'].centroid_values_gpu, 
                        self.quantities['ymomentum'].centroid_values_gpu,
                        block = (W1, W2, W3),
                        grid=((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                        )
        else:
            Domain.protect_against_infinitesimal_and_negative_heights(self)

    def extrapolate_second_order_sw(self):
        if  self.using_gpu:
            W1 = 8
            W2 = 1
            W3 = 1

            if self.extrapolate_velocity_second_order == 1:
                self.extrapolate_second_order_sw_true_func(
                    numpy.float64(self.epsilon),
                    numpy.float64(self.minimum_allowed_height),
                    numpy.float64(self.beta_w),
                    numpy.float64(self.beta_w_dry),
                    numpy.float64(self.beta_uh),
                    numpy.float64(self.beta_uh_dry),
                    numpy.float64(self.beta_vh),
                    numpy.float64(self.beta_vh_dry),
                    numpy.float64(self.optimise_dry_cells),

                    self.surrogate_neighbours_gpu,
                    self.number_of_boundaries_gpu,
                    self.centroid_coordinates_gpu,
                    self.quantities['stage'].centroid_values_gpu,
                    self.quantities['elevation'].centroid_values_gpu,
                    self.quantities['xmomentum'].centroid_values_gpu,
                    self.quantities['ymomentum'].centroid_values_gpu,
                    self.vertex_coordinates_gpu,
                    self.quantities['stage'].vertex_values_gpu,
                    self.quantities['elevation'].vertex_values_gpu,
                    self.quantities['xmomentum'].vertex_values_gpu,
                    self.quantities['ymomentum'].vertex_values_gpu,
                    self.stage_centroid_store_gpu,
                    self.xmomentum_centroid_store_gpu,
                    self.ymomentum_centroid_store_gpu,
                    block = (W1, W2, W3),
                    grid=((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                    )
            else:
                self.extrapolate_second_order_sw_false_func(
                    numpy.float64(self.epsilon),
                    numpy.float64(self.minimum_allowed_height),
                    numpy.float64(self.beta_w),
                    numpy.float64(self.beta_w_dry),
                    numpy.float64(self.beta_uh),
                    numpy.float64(self.beta_uh_dry),
                    numpy.float64(self.beta_vh),
                    numpy.float64(self.beta_vh_dry),
                    numpy.float64(self.optimise_dry_cells),
                    self.surrogate_neighbours_gpu,
                    self.number_of_boundaries_gpu,
                    self.centroid_coordinates_gpu,
                    self.quantities['stage'].centroid_values,
                    self.quantities['elevation'].centroid_values,
                    self.quantities['xmomentum'].centroid_values,
                    self.quantities['ymomentum'].centroid_values,
                    self.vertex_coordinates,
                    self.quantities['stage'].vertex_values,
                    self.quantities['elevation'].vertex_values,
                    self.quantities['xmomentum'].vertex_values,
                    self.quantities['ymomentum'].vertex_values,
                    self.stage_centroid_store_gpu,
                    self.xmomentum_centroid_store_gpu,
                    self.ymomentum_centroid_store_gpu,
                    block = (W1, W2, W3),
                    grid=((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3),1)
                    )
        else:
            Domain.extrapolate_second_order_sw(self)

    def evolve(self, 
                yieldstep=None,
                finaltime=None,
                duration=None,
                skip_initial_step=False):

        if self.using_gpu:
            #self.lock_array_page()

            self.allocate_device_array()

            self.timestep_array[:] = self.evolve_max_timestep
            
            self.asynchronous_transfer()



        for t in Domain.evolve(self, yieldstep=yieldstep,
                    finaltime=finaltime,
                    duration=duration,
                    skip_initial_step=skip_initial_step):

            yield(t)


        global auto_init_context
        if not auto_init_context:
            ctx.pop()






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
    assert numpy.allclose(a, temp_page_lock_p)           
    return temp_page_lock_p

def get_device_array(a):
    return drv.mem_alloc(a.nbytes)

def asy_cpy(a, a_gpu):
    strm = drv.Stream()
    drv.memcpy_htod_async(a_gpu, a, strm)
    return strm
