#!/usr/bin/env python

# System module
import numpy 
import sys
import time



# ANUGA module
from anuga import Domain
from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import \
        Transmissive_boundary, \
        Dirichlet_boundary, \
        Compute_fluxes_boundary, \
        Time_boundary, Time_boundary, \
        Time_boundary



# PyCUDA module
import pycuda.driver as drv
from pycuda.compiler import SourceModule
auto_init_context = True
using_page_locked = False
if auto_init_context:
    import pycuda.autoinit
    dev = pycuda.autoinit.device
    ctx = pycuda.autoinit.context
else:
    drv.init()
    dev = drv.Device(0)
    ctx = dev.make_context( drv.ctx_flags.MAP_HOST)



# Config 
#from anuga_cuda import kernel_path as kp
#from anuga_cuda import kernel_block_configuration as kbc
from anuga_cuda import *
kp = kernel_path
kbc = kernel_block_configuration 



class GPU_domain(Domain):
    def __init__(self, 
            coordinates=None,
            vertices=None,
            boundary=None,
            source=None,
            triangles=None,
            conserved_quantities=None,
            evolved_quantities=None,
            other_quantities=None,
            tagged_elements=None,
            geo_reference=None,
            use_inscribed_circle=False,
            mesh_filename=None,
            use_cache=False,
            verbose=False,
            full_send_dict=None,
            ghost_recv_dict=None,
            starttime=0.0,
            processor=0,
            numproc=1,
            number_of_full_nodes=None,
            number_of_full_triangles=None,
            ghost_layer_width=2,
            using_gpu=False,
            cotesting=False,
            stream= False,
            rearranged=False,
            domain=None): 

        if domain == None:
            Domain.__init__(self,
                coordinates,
                vertices,
                boundary,
                tagged_elements,
                geo_reference,
                use_inscribed_circle,
                mesh_filename,
                use_cache,
                verbose,
                conserved_quantities ,
                evolved_quantities ,
                other_quantities ,
                full_send_dict,
                ghost_recv_dict,
                starttime,
                processor,
                numproc,
                number_of_full_nodes,
                number_of_full_triangles,
                ghost_layer_width)
        else:
            self.__dict__.update(domain.__dict__)

        self.boundary_index = {}


        self.using_gpu = using_gpu
        self.cotesting = cotesting
        self.rearranged_domain = rearranged



        print '\n --> Device attributes'
        print '    Name:', dev.name()
        print '    Compute capability:', dev.compute_capability()
        print '    Concurrent Kernels:', \
            bool(dev.get_attribute(
                drv.device_attribute.CONCURRENT_KERNELS))
        print "    Current cache/shared memory configure is ", \
            ctx.get_cache_config()
        ctx.set_cache_config(drv.func_cache.PREFER_L1)
        #ctx.set_cache_config(drv.func_cache.PREFER_SHARED)
        #ctx.set_cache_config(drv.func_cache.PREFER_EQUAL)
        #ctx.set_cache_config(drv.func_cache.PREFER_NONE)


        if stream:
            if bool(dev.get_attribute(
                drv.device_attribute.CONCURRENT_KERNELS)):
                self.using_stream = stream
            else:
                print ' *** Disable stream ***'
                self.using_stream = False
        else:
            self.using_stream = False

        #self.end_event = drv.Event()


        #if not cotesting:
        #    self.equip_kernel_functions()

        if using_page_locked:
            self.timestep_array = drv.pagelocked_zeros(
                    self.number_of_elements,
                    dtype = numpy.float64,
                    mem_flags=drv.host_alloc_flags.DEVICEMAP)
        else:
            self.timestep_array = numpy.zeros(
                    self.number_of_elements,
                    dtype = numpy.float64)
 

    """ GPU_domain functions
    """


    def equip_kernel_functions(self):
        """ Get kernel functions
        """

        # compute_fluxes function
        self.compute_fluxes_mod = get_sourceModule( 
                kp["compute_fluxes_dir"],
                "compute_fluxes.cu", self.rearranged_domain)

        self.compute_fluxes_func = self.compute_fluxes_mod.get_function(
                #"compute_fluxes_central_structure_cuda_single")
                "compute_fluxes_central_structure_CUDA")
        #self.compute_fluxes_func.set_cache_config(
        #        drv.func_cache.PREFER_L1)
        self.compute_fluxes_func.prepare("idddddiPPPPPPPPPPPPPPPPPPP")
        self.compute_fluxes_central_structure_block = \
                kbc["compute_fluxes_fun"]
        


        # gravity_wb function
        self.gravity_wb_mod = get_sourceModule( kp["gravity_dir"], 
                "gravity.cu", self.rearranged_domain)

        self.gravity_wb_func = self.gravity_wb_mod.get_function(
                "gravity_wb")
        self.gravity_wb_func.prepare("idPPPPPPPPPPP")
        self.gravity_wb_block = kbc["gravity_fun"]



        # extrapolate_second_order_sw function
        self.extrapolate_second_order_sw_mod = get_sourceModule(
                kp["extrapolate_dir"], "extrapolate_second_order_sw.cu", 
                self.rearranged_domain)

        self.extrapolate_second_order_sw_true_func = \
                self.extrapolate_second_order_sw_mod.get_function(
                        "extrapolate_second_order_sw_true")
        self.extrapolate_second_order_sw_true_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.extrapolate_second_order_sw_true_func.prepare(
                "iddddddddiPPPPPPPPPPPP")

        self.extrapolate_second_order_sw_false_func = \
                self.extrapolate_second_order_sw_mod.get_function(
                        "extrapolate_second_order_sw_false")
        self.extrapolate_second_order_sw_false_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.extrapolate_second_order_sw_false_func.prepare(
                "iddddddddiPPPPPPPPPPPP")

        self.extrapolate_velocity_second_order_true_func = \
                self.extrapolate_second_order_sw_mod.get_function(
                        "extrapolate_velocity_second_order_true")
        self.extrapolate_velocity_second_order_true_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.extrapolate_velocity_second_order_true_func.prepare(
                "idPPPPPP")


        self.extrapolate_second_order_sw_true_block = \
            kbc["extrapolate_second_order_sw_fun"]
        self.extrapolate_second_order_sw_false_block= \
            kbc["extrapolate_second_order_sw_fun"]
        self.extrapolate_velocity_second_order_true_block = \
            kbc["extrapolate_velocity_second_order_true_fun"]



        # extrapolate_second_order_and_limit_by_vertex_or_edge function
        self.extrapolate_second_order_and_limit_by_vertex_or_edge_mod = \
                get_sourceModule( kp["extrapolate_dir"],
                "extrapolate_second_order_and_limit_by_vertex_or_edge.cu",
                self.rearranged_domain)

        self.extrapolate_second_order_and_limit_by_vertex_func = \
            self.extrapolate_second_order_and_limit_by_vertex_or_edge_mod.\
                    get_function(
                    "extrapolate_second_order_and_limit_by_vertex")
        self.extrapolate_second_order_and_limit_by_vertex_func.\
                set_cache_config(drv.func_cache.PREFER_L1)
        self.extrapolate_second_order_and_limit_by_vertex_func.prepare(
                "idPPPPPPPPPP")


        self.extrapolate_second_order_and_limit_by_edge_func = \
            self.extrapolate_second_order_and_limit_by_vertex_or_edge_mod.\
                get_function(
                    "extrapolate_second_order_and_limit_by_edge")
        self.extrapolate_second_order_and_limit_by_edge_func.\
                set_cache_config(drv.func_cache.PREFER_L1)
        self.extrapolate_second_order_and_limit_by_edge_func.prepare(
                "idPPPPPPPPPPP")

        self.extrapolate_second_order_and_limit_by_vertex_block = \
            kbc["extrapolate_second_order_and_limit_by_vertex_fun"]
        self.extrapolate_second_order_and_limit_by_edge_block = \
            kbc["extrapolate_second_order_and_limit_by_edge_fun"]



        # extrapolate_first_order function
        self.extrapolate_first_order_mod = get_sourceModule(
                kp["extrapolate_dir"], "extrapolate_first_order.cu",
                self.rearranged_domain)

        self.extrapolate_first_order_func = \
                self.extrapolate_first_order_mod.get_function(
                        "extrapolate_first_order")
        self.extrapolate_first_order_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.extrapolate_first_order_func.prepare( "iPPP" )
        self.extrapolate_first_order_block = \
                kbc["extrapolate_first_order_fun"]



        # extrapolate_second_order_edge_swb2 function
        self.extrapolate_second_order_edge_swb2_mod = get_sourceModule(
                kp["extrapolate_dir"],
                "extrapolate_second_order_edge_swb2.cu", 
                self.rearranged_domain)

        self.extrapolate_second_order_edge_swb2_func = \
            self.extrapolate_second_order_edge_swb2_mod.get_function(
                    "_extrapolate_second_order_edge_sw")
        self.extrapolate_second_order_edge_swb2_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.extrapolate_second_order_edge_swb2_func.prepare(
                "iiiddddddddPPPPPPPPPPPPPPPPPPPPPP")

        self.extrapolate_second_order_edge_swb2_block =\
            kbc["extrapolate_second_order_edge_swb2_fun"]



        # protect function
        self.protect_mod = get_sourceModule(kp["protect_dir"], 
                "protect.cu", self.rearranged_domain)

        self.protect_sw_func = self.protect_mod.get_function("_protect_sw")
        self.protect_sw_func.set_cache_config(drv.func_cache.PREFER_L1)
        self.protect_sw_func.prepare( "idddPPPP" )

        self.protect_swb2_func=self.protect_mod.get_function(
                "_protect_swb2")
        self.protect_swb2_func.set_cache_config(drv.func_cache.PREFER_L1)
        self.protect_swb2_func.prepare( "idddPPPPPPP" )

        self.protect_sw_block = kbc["protect_sw_ext_fun"]
        self.protect_swb2_block= kbc["protect_swb2_fun"]



        # balance function
        self.balance_deep_and_shallow_mod = get_sourceModule(
                kp["balance_dir"], "balance_deep_and_shallow.cu", 
                self.rearranged_domain)

        self.balance_deep_and_shallow_func = \
                self.balance_deep_and_shallow_mod.get_function(
                        "_balance_deep_and_shallow")
        self.balance_deep_and_shallow_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.balance_deep_and_shallow_func.prepare( "iddiiPPPPPPPP" )

        self.balance_deep_and_shallow_block = kbc["balance_fun"]
            


        # interpolate_from_vertices_to_edges function
        self.interpolate_from_vertices_to_edges_mod = get_sourceModule(
                kp["interpolate_dir"],
                "interpolate_from_vertices_to_edges.cu", 
                self.rearranged_domain)

        self.interpolate_from_vertices_to_edges_func = \
                self.interpolate_from_vertices_to_edges_mod.get_function(
                        "_interpolate_from_vertices_to_edges")
        self.interpolate_from_vertices_to_edges_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.interpolate_from_vertices_to_edges_func.prepare( "iPP" )

        self.interpolate_from_vertices_to_edges_block = \
                kbc["interpolate_fun"]



        # evaluate_segment function
        self.evaluate_segment_mod = get_sourceModule(
                kp["evaluate_dir"], "evaluate_segment.cu", 
                self.rearranged_domain)

        self.evaluate_segment_reflective_func = \
                self.evaluate_segment_mod.get_function(
                        "evaluate_segment_reflective")
        self.evaluate_segment_reflective_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.evaluate_segment_reflective_func.prepare( 
                "iiPPPPPPPPPPPPPPPPPP")

        self.evaluate_segment_dirichlet_1_func = \
                self.evaluate_segment_mod.get_function(
                        "evaluate_segment_dirichlet_1")
        self.evaluate_segment_dirichlet_1_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.evaluate_segment_dirichlet_1_func.prepare( "iPPPPP " )

        self.evaluate_segment_dirichlet_2_func = \
                self.evaluate_segment_mod.get_function(
                        "evaluate_segment_dirichlet_2")
        self.evaluate_segment_dirichlet_2_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.evaluate_segment_dirichlet_2_func.prepare( "idPP" )

        self.evaluate_segment_reflective_block = \
            kbc["evaluate_segment_reflective_fun"]
        self.evaluate_segment_dirichlet_1_block = \
            kbc["evaluate_segment_dirichlet_1_fun"]
        self.evaluate_segment_dirichlet_2_block = \
            kbc["evaluate_segment_dirichlet_2_fun"]
            


        # get_absolute function
        self.get_absolute_mod = get_sourceModule(
                kp["get_absolute_dir"], "get_absolute.cu", 
                self.rearranged_domain)
            
        self.get_absolute_func = \
                self.get_absolute_mod.get_function("get_absolute")
        self.get_absolute_func.set_cache_config(drv.func_cache.PREFER_L1)
        self.get_absolute_func.prepare( "iddP" )

        self.get_absolute_block = kbc["get_absolute_fun"]



        # manning_friction function
        self.manning_friction_mod = get_sourceModule(
                kp["manning_friction_dir"], "manning_friction.cu", 
                self.rearranged_domain)

        self.manning_friction_sloped_func = \
                self.manning_friction_mod.get_function(
                        "_manning_friction_sloped")
        self.manning_friction_sloped_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.manning_friction_sloped_func.prepare( "iddPPPPPPPP" )

        self.manning_friction_flat_func = \
                self.manning_friction_mod.get_function(
                        "_manning_friction_flat")
        self.manning_friction_flat_func.set_cache_config(
                drv.func_cache.PREFER_L1)
        self.manning_friction_flat_func.prepare( "iddPPPPPPP" ) 

        self.manning_friction_sloped_blcok = \
                kbc["manning_friction_sloped_fun"]
        self.manning_friction_flat_block = \
                kbc["manning_friction_flat_fun"]



        # saxpy_centroid_values function
        self.saxpy_centroid_values_mod = get_sourceModule(
                kp["saxpy_dir"], "saxpy_centroid_values.cu", 
                self.rearranged_domain)

        self.saxpy_centroid_values_func = \
                self.saxpy_centroid_values_mod.get_function(
                        "_saxpy_centroid_values")
        self.saxpy_centroid_values_func.prepare( "iddPP" )

        self.saxpy_centroid_values_block = kbc["saxpy_fun"]



        # set_boundry function
        self.set_boundary_mod = get_sourceModule( kp["set_boundary_dir"], 
                "set_boundary.cu", self.rearranged_domain)
                
        self.set_boundary_values_from_edges_func = \
                self.set_boundary_mod.get_function(
                    "set_boundary_values_from_edges")
        self.set_boundary_values_from_edges_func.prepare( "iPPPP" )

        self.set_boundary_block = \
            kbc["set_boundary_values_from_edges_fun"]



        # update_centroids_of_velocities_and_height function
        self.update_centroids_of_velocities_and_height_mod = \
                get_sourceModule( kp["update_centroids_dir"],
                        "update_centroids_of_velocities_and_height.cu", 
                        self.rearranged_domain)
        
        self.update_centroids_of_velocities_and_height_func = \
            self.update_centroids_of_velocities_and_height_mod.\
            get_function("update_centroids_of_velocities_and_height")
        self.update_centroids_of_velocities_and_height_func.prepare(
                "iiPPPPPPPPPPPPPP")
        
        self.update_centroids_of_velocities_and_height_block = \
                kbc["update_centroids_fun"]
            

        
        # update function
        self.update_mod = get_sourceModule( kp["update_dir"], 
                "update.cu", self.rearranged_domain)
        
        self.update_func = self.update_mod.get_function("update")
        self.update_func.prepare( "idPPP" )

        self.update_block = kbc["update_fun"]



    def lock_host_page(self):
        """ Use page-locked memory

        Register host pageable memory to lock their page.
        """

        self.neighbours = get_page_locked_array(self.neighbours)
        self.neighbour_edges = get_page_locked_array(self.neighbour_edges)
        self.normals = get_page_locked_array(self.normals)
        self.edgelengths = get_page_locked_array(self.edgelengths)
        self.radii = get_page_locked_array(self.radii)
        self.areas = get_page_locked_array(self.areas)
        self.tri_full_flag = get_page_locked_array(self.tri_full_flag)
        self.max_speed = get_page_locked_array(self.max_speed)
        self.vertex_coordinates = \
            get_page_locked_array(self.vertex_coordinates)

        
        
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
        """ Allocate device memory and copy data to device
        """

        # auxiliary arrays
        self.timestep_array_gpu = get_device_array(self.timestep_array)
        self.stage_centroid_store_gpu = get_device_array(
                self.quantities['stage'].centroid_values)
        self.xmomentum_centroid_store_gpu = get_device_array(
                self.quantities['xmomentum'].centroid_values)
        self.ymomentum_centroid_store_gpu = get_device_array(
                self.quantities['ymomentum'].centroid_values)
        self.mim_elevation_edgevalue_gpu = get_device_array(
                self.quantities['ymomentum'].centroid_values)
        self.max_elevation_edgevalue_gpu = get_device_array(
                self.quantities['ymomentum'].centroid_values)
        self.count_wet_neighbours_gpu =  drv.mem_alloc(
                self.number_of_elements * numpy.int32().nbytes)
        
        for tag in self.tag_boundary_cells:
            if self.boundary_map[tag] is None:
                continue
            self.boundary_index[tag] = []

            ids = self.tag_boundary_cells[tag] 
            
            self.boundary_index[tag].append( 
                    get_device_array( numpy.asarray(ids) ))
            self.boundary_index[tag].append(
                    get_device_array( numpy.asarray(self.boundary_cells[ids])) )
            self.boundary_index[tag].append(
                    get_device_array( numpy.asarray(self.boundary_edges[ids])) )


        # domain arrays
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
        self.boundary_cells_gpu = get_device_array(
                numpy.asarray(self.boundary_cells))
        self.boundary_edges_gpu = get_device_array(
                numpy.asarray(self.boundary_edges))

        
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
        self.quantities['height'].edge_values_gpu = \
            get_device_array(self.quantities['height'].edge_values)
        self.quantities['xvelocity'].edge_values_gpu = \
            get_device_array(self.quantities['xvelocity'].edge_values)
        self.quantities['yvelocity'].edge_values_gpu = \
            get_device_array(self.quantities['yvelocity'].edge_values)
        
        # get boundary values
        self.quantities['stage'].boundary_values_gpu = \
            get_device_array(self.quantities['stage'].boundary_values)
        self.quantities['elevation'].boundary_values_gpu = \
            get_device_array(self.quantities['elevation'].boundary_values)
        self.quantities['height'].boundary_values_gpu = \
            get_device_array(self.quantities['height'].boundary_values)
        self.quantities['xmomentum'].boundary_values_gpu = \
            get_device_array( self.quantities['xmomentum'].boundary_values)
        self.quantities['ymomentum'].boundary_values_gpu = \
            get_device_array(self.quantities['ymomentum'].boundary_values)
        self.quantities['xvelocity'].boundary_values_gpu = \
            get_device_array( self.quantities['xvelocity'].boundary_values)
        self.quantities['yvelocity'].boundary_values_gpu = \
            get_device_array(self.quantities['yvelocity'].boundary_values)
        
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
            get_device_array(self.quantities['elevation'].vertex_values)
        self.quantities['xmomentum'].vertex_values_gpu = \
            get_device_array(self.quantities['xmomentum'].vertex_values)
        self.quantities['ymomentum'].vertex_values_gpu = \
            get_device_array(self.quantities['ymomentum'].vertex_values)
        self.quantities['height'].vertex_values_gpu = \
            get_device_array(self.quantities['height'].vertex_values)
        self.quantities['xvelocity'].vertex_values_gpu = \
            get_device_array(self.quantities['xvelocity'].vertex_values)
        self.quantities['yvelocity'].vertex_values_gpu = \
            get_device_array(self.quantities['yvelocity'].vertex_values)
        
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
        self.quantities['friction'].centroid_values_gpu = \
            get_device_array(self.quantities['friction'].centroid_values)
        self.quantities['height'].centroid_values_gpu = \
            get_device_array(self.quantities['height'].centroid_values)
        self.quantities['xvelocity'].centroid_values_gpu = \
            get_device_array(self.quantities['xvelocity'].centroid_values)
        self.quantities['yvelocity'].centroid_values_gpu = \
            get_device_array(self.quantities['yvelocity'].centroid_values)
        
        for name in self.conserved_quantities:
            """ ['stage', 'xmomentum', 'ymomentum']"""
            Q = self.quantities[name]
            Q.x_gradient_gpu = get_device_array(Q.x_gradient)
            Q.y_gradient_gpu = get_device_array(Q.y_gradient)
            Q.centroid_backup_values_gpu = \
                    get_device_array(Q.centroid_values)
            Q.semi_implicit_update_gpu = \
                    get_device_array( Q.semi_implicit_update)
        
        for name in ['height', 'xvelocity', 'yvelocity']:
            Q = self.quantities[name]
            Q.x_gradient_gpu = get_device_array(Q.x_gradient)
            Q.y_gradient_gpu = get_device_array(Q.y_gradient)



    def copy_back_necessary_data(self):
        pass       


    def asynchronous_transfer(self):
        """ Asynchronous transfer from host to device

        Asynchronous transfer data from page-locked memory can hide 
        transmission overheads by proceeding several transfer together and 
        overlapping with kernel executing time
        """

        # auxiliary arrays
        for tag in self.tag_boundary_cells:
            if self.boundary_map is None:
                continue
            ids = self.tag_boundary_cells[tag] 
            asy_cpy( numpy.asarray(ids), self.boundary_index[tag][0])
            asy_cpy( numpy.asarray(self.boundary_cells[ids]), 
                    self.boundary_index[tag][1])
            asy_cpy( numpy.asarray(self.boundary_edges[ids]),
                    self.boundary_index[tag][2])

        asy_cpy(self.timestep_array, self.timestep_array_gpu)

        # domain arrays
        asy_cpy(self.neighbours, self.neighbours_gpu)
        asy_cpy(self.neighbour_edges, self.neighbour_edges_gpu)
        asy_cpy(self.surrogate_neighbours, self.surrogate_neighbours_gpu)
        asy_cpy(self.normals, self.normals_gpu)
        asy_cpy(self.edgelengths, self.edgelengths_gpu)
        asy_cpy(self.radii, self.radii_gpu)
        asy_cpy(self.areas, self.areas_gpu)
        asy_cpy(self.tri_full_flag, self.tri_full_flag_gpu)
        asy_cpy(self.max_speed, self.max_speed_gpu)
        asy_cpy(self.boundary_cells, self.boundary_cells_gpu)
        asy_cpy(self.boundary_edges, self.boundary_edges_gpu)

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
        asy_cpy(
                self.quantities['height'].edge_values,
                self.quantities['height'].edge_values_gpu)
        asy_cpy(
                self.quantities['xvelocity'].edge_values,
                self.quantities['xvelocity'].edge_values_gpu)
        asy_cpy(
                self.quantities['yvelocity'].edge_values,
                self.quantities['yvelocity'].edge_values_gpu)

        # get boundary values
        asy_cpy(
                self.quantities['stage'].boundary_values,
                self.quantities['stage'].boundary_values_gpu)
        asy_cpy(
                self.quantities['elevation'].boundary_values,
                self.quantities['elevation'].boundary_values_gpu)
        asy_cpy(
                self.quantities['height'].boundary_values,
                self.quantities['height'].boundary_values_gpu)
        asy_cpy(
                self.quantities['xmomentum'].boundary_values,
                self.quantities['xmomentum'].boundary_values_gpu)
        asy_cpy(
                self.quantities['ymomentum'].boundary_values,
                self.quantities['ymomentum'].boundary_values_gpu)
        asy_cpy(
                self.quantities['xvelocity'].boundary_values,
                self.quantities['xvelocity'].boundary_values_gpu)
        asy_cpy(
                self.quantities['yvelocity'].boundary_values,
                self.quantities['yvelocity'].boundary_values_gpu)

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

        # semi_implicit_update
        asy_cpy(
                self.quantities['stage'].semi_implicit_update,
                self.quantities['stage'].semi_implicit_update_gpu)
        asy_cpy(
                self.quantities['xmomentum'].semi_implicit_update,
                self.quantities['xmomentum'].semi_implicit_update_gpu)
        asy_cpy(
                self.quantities['ymomentum'].semi_implicit_update,
                self.quantities['ymomentum'].semi_implicit_update_gpu)


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
        asy_cpy(
                self.quantities['friction'].centroid_values,
                self.quantities['friction'].centroid_values_gpu)
        asy_cpy(
                self.quantities['height'].centroid_values,
                self.quantities['height'].centroid_values_gpu)
        asy_cpy(
                self.quantities['xvelocity'].centroid_values,
                self.quantities['xvelocity'].centroid_values_gpu)
        asy_cpy(
                self.quantities['yvelocity'].centroid_values,
                self.quantities['yvelocity'].centroid_values_gpu)


        # conserved_quantities arrays
        for name in self.conserved_quantities:
            Q = self.quantities[name]
            asy_cpy( Q.x_gradient, Q.x_gradient_gpu)
            asy_cpy( Q.y_gradient, Q.y_gradient_gpu)

        for name in ['height', 'xvelocity', 'yvelocity']:
            Q = self.quantities[name]
            asy_cpy( Q.vertex_values, Q.vertex_values_gpu)
            asy_cpy( Q.x_gradient, Q.x_gradient_gpu)
            asy_cpy( Q.x_gradient, Q.x_gradient_gpu)
            asy_cpy( Q.y_gradient, Q.y_gradient_gpu)





    """ Reloading functions from ANUGA
    """



    # 5th level cotesting
    def apply_protection_against_isolated_degenerate_timesteps(self):
        if self.using_gpu:
            if self.protect_against_isolated_degenerate_timesteps is False:
                return
            
            cpy_back(self.max_speed, self.max_speed_gpu)
            ctx.synchronize()

            if numpy.max(self.max_speed) < 10.0:
                return

            from anuga.utilities.numerical_tools import \
                    histogram, create_bins
            
            bins = create_bins(self.max_speed, 10)
            hist = histogram(self.max_speed, bins)

            if len(hist) > 1 and hist[-1] > 0 and \
                hist[4] == hist[5] == hist[6] == hist[7] == hist[8] == 0:

                d = 0
                for i in range(self.number_of_triangles):
                    if self.max_speed[i] > bins[-1]:
                        print "1"
                        self.get_quantity('xmomentum').\
                            set_values(0.0, indices=[i])
                        self.get_quantity('ymomentum').\
                            set_values(0.0, indices=[i])
                        self.max_speed[i]=0.0
                        d += 1
        else:
            Domain.apply_protection_against_isolated_degenerate_timesteps(self)



    # 4th level cotesting
    # Using Stream
    def balance_deep_and_shallow(self):
        if  self.using_gpu:
            N = self.number_of_elements
            W1 = self.balance_deep_and_shallow_block
            W2 = 1
            W3 = 1

            s = self.quantities['stage']
            b = self.quantities['elevation']
            x = self.quantities['xmomentum']
            y = self.quantities['ymomentum']

            s.centroid_values_gpu=get_device_array(s.centroid_values, True)
            b.centroid_values_gpu=get_device_array(b.centroid_values, True)
            x.centroid_values_gpu=get_device_array(x.centroid_values, True)
            y.centroid_values_gpu=get_device_array(y.centroid_values, True)
            

            s.vertex_values_gpu =get_device_array(s.vertex_values, True)
            b.vertex_values_gpu =get_device_array(b.vertex_values, True)
            x.vertex_values_gpu =get_device_array(x.vertex_values, True)
            y.vertex_values_gpu =get_device_array(y.vertex_values, True)
            

            self.balance_deep_and_shallow_func.prepared_call(
                ((N+W1*W2*W3-1)/(W1*W2*W3), 1),
                (W1, W2, W3),
                #self.stream[balance_stream],
                N,
                self.H0,
                self.alpha_balance,

                self.tight_slope_limiters,
                self.use_centroid_velocities,
                
                self.quantities['stage'].centroid_values_gpu,
                self.quantities['elevation'].centroid_values_gpu,
                self.quantities['stage'].vertex_values_gpu,
                self.quantities['elevation'].vertex_values_gpu,

                self.quantities['xmomentum'].centroid_values_gpu,
                self.quantities['ymomentum'].centroid_values_gpu,
                self.quantities['xmomentum'].vertex_values_gpu,
                self.quantities['ymomentum'].vertex_values_gpu
                )
            
            cpy_back(s.centroid_values, s.centroid_values_gpu, False)
            cpy_back(b.centroid_values, b.centroid_values_gpu, False)
            cpy_back(x.centroid_values, x.centroid_values_gpu, False)
            cpy_back(y.centroid_values, y.centroid_values_gpu, False)
            
            cpy_back(s.vertex_values, s.vertex_values_gpu, False)
            cpy_back(b.vertex_values, b.vertex_values_gpu, False)
            cpy_back(x.vertex_values, x.vertex_values_gpu, False)
            cpy_back(y.vertex_values, y.vertex_values_gpu, False)
            
            s.centroid_values_gpu.free()
            b.centroid_values_gpu.free()
            x.centroid_values_gpu.free()
            y.centroid_values_gpu.free()
            
            s.vertex_values_gpu.free() 
            b.vertex_values_gpu.free()
            x.vertex_values_gpu.free()
            y.vertex_values_gpu.free()
        else:
            Domain.balance_deep_and_shallow(self)
        if self.cotesting:
            Domain.balance_deep_and_shallow(self.cotesting_domain)
            test_balance_deep_and_shallow(self)





    # 4th level cotesting
    # Using Stream
    def protect_against_infinitesimal_and_negative_heights(self):
        if  self.using_gpu:
            N = self.number_of_elements
            W2 = 1
            W3 = 1

            if self.flow_algorithm == 'tsunami':
                W1 = self.protect_swb2_fun
                self.protect_swb2_func(
                        numpy.int32(N),
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
                        grid=((N+W1*W2*W3-1)/(W1*W2*W3),1),
                        stream = self.stream[ protect_swb2_stream ]
                        )
                    
            else:
                s = self.quantities['stage']
                b = self.quantities['elevation']
                x = self.quantities['xmomentum']
                y = self.quantities['ymomentum']

                s.centroid_values_gpu=get_device_array(s.centroid_values, True)
                b.centroid_values_gpu=get_device_array(b.centroid_values, True)
                x.centroid_values_gpu=get_device_array(x.centroid_values, True)
                y.centroid_values_gpu=get_device_array(y.centroid_values, True)

                W1 = self.protect_sw_block
                self.protect_sw_func.prepared_call(
                        ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                        (W1, W2, W3),
                        #self.stream[ protect_sw_stream ],
                        N,
                        self.minimum_allowed_height, 
                        self.maximum_allowed_speed,
                        self.epsilon, 
                        self.quantities['stage'].centroid_values_gpu,
                        self.quantities['elevation'].centroid_values_gpu, 
                        self.quantities['xmomentum'].centroid_values_gpu, 
                        self.quantities['ymomentum'].centroid_values_gpu
                        )

                cpy_back(s.centroid_values, s.centroid_values_gpu, False)
                cpy_back(b.centroid_values, b.centroid_values_gpu, False)
                cpy_back(x.centroid_values, x.centroid_values_gpu, False)
                cpy_back(y.centroid_values, y.centroid_values_gpu, False)

                s.centroid_values_gpu.free()
                b.centroid_values_gpu.free()
                x.centroid_values_gpu.free()
                y.centroid_values_gpu.free()
        else:
            Domain.protect_against_infinitesimal_and_negative_heights(self)

        if self.cotesting:
            Domain.protect_against_infinitesimal_and_negative_heights(
                    self.cotesting_domain)
            test_protect_against_infinitesimal_and_negative_heights(self)
            

    # 4th level cotesting
    # Using Stream
    def extrapolate_second_order_sw(self):
        if  self.using_gpu:
            N = self.number_of_elements
            W2 = 1
            W3 = 1

            if self.extrapolate_velocity_second_order :
                s = self.quantities['stage']
                b = self.quantities['elevation']
                x = self.quantities['xmomentum']
                y = self.quantities['ymomentum']

                self.surrogate_neighbours_gpu=get_device_array(self.surrogate_neighbours,True)
                self.number_of_boundaries_gpu=get_device_array(self.number_of_boundaries,True)
                self.centroid_coordinates_gpu=get_device_array(self.centroid_coordinates,True)
                self.vertex_coordinates_gpu=get_device_array(self.vertex_coordinates,True)
                
                s.centroid_values_gpu=get_device_array(s.centroid_values, True)
                b.centroid_values_gpu=get_device_array(b.centroid_values, True)
                x.centroid_values_gpu=get_device_array(x.centroid_values, True)
                y.centroid_values_gpu=get_device_array(y.centroid_values, True)
                s.vertex_values_gpu=get_device_array(s.vertex_values, True)
                b.vertex_values_gpu=get_device_array(b.vertex_values, True)
                x.vertex_values_gpu=get_device_array(x.vertex_values, True)
                y.vertex_values_gpu=get_device_array(y.vertex_values, True)
                self.xmomentum_centroid_store_gpu=get_device_array(x.centroid_values, False)
                self.ymomentum_centroid_store_gpu=get_device_array(x.centroid_values, False)

                W1 = self.extrapolate_velocity_second_order_true_block
                self.extrapolate_velocity_second_order_true_func.\
                    prepared_call(
                        ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                        (W1, W2, W3),
                        #self.stream[ extra_2_sw_stream]
                        N,
                        self.minimum_allowed_height,

                        self.quantities['stage'].centroid_values_gpu,
                        self.quantities['elevation'].centroid_values_gpu,
                        self.quantities['xmomentum'].centroid_values_gpu,
                        self.xmomentum_centroid_store_gpu,
                        self.quantities['ymomentum'].centroid_values_gpu,
                        self.ymomentum_centroid_store_gpu
                        )
                
                W1 = self.extrapolate_second_order_sw_true_block
                self.extrapolate_second_order_sw_true_func.prepared_call(
                    ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                    (W1, W2, W3),
                    #self.stream[ extra_2_sw_stream],
                    N,
                    self.epsilon,
                    self.minimum_allowed_height,
                    self.beta_w,
                    self.beta_w_dry,
                    self.beta_uh,
                    self.beta_uh_dry,
                    self.beta_vh,
                    self.beta_vh_dry,
                    self.optimise_dry_cells,

                    self.surrogate_neighbours_gpu,
                    self.number_of_boundaries_gpu,
                    self.centroid_coordinates_gpu,
                    self.quantities['stage'].centroid_values_gpu,
                    self.quantities['elevation'].centroid_values_gpu,
                    #self.quantities['xmomentum'].centroid_values_gpu,
                    #self.quantities['ymomentum'].centroid_values_gpu,
                    self.xmomentum_centroid_store_gpu,
                    self.ymomentum_centroid_store_gpu,
                    self.vertex_coordinates_gpu,
                    self.quantities['stage'].vertex_values_gpu,
                    self.quantities['elevation'].vertex_values_gpu,
                    self.quantities['xmomentum'].vertex_values_gpu,
                    self.quantities['ymomentum'].vertex_values_gpu
                    )

                self.surrogate_neighbours_gpu.free()
                self.number_of_boundaries_gpu.free()
                self.centroid_coordinates_gpu.free()
                self.vertex_coordinates_gpu.free()

                self.xmomentum_centroid_store_gpu.free()
                self.ymomentum_centroid_store_gpu.free()

                cpy_back(s.vertex_values, s.vertex_values_gpu, False)
                cpy_back(b.vertex_values, b.vertex_values_gpu, False)
                cpy_back(x.vertex_values, x.vertex_values_gpu, False)
                cpy_back(y.vertex_values, y.vertex_values_gpu, False)

                s.vertex_values_gpu.free()
                b.vertex_values_gpu.free()
                x.vertex_values_gpu.free()
                y.vertex_values_gpu.free()

                cpy_back(s.centroid_values, s.centroid_values_gpu, False)
                cpy_back(b.centroid_values, b.centroid_values_gpu, False)
                cpy_back(x.centroid_values, x.centroid_values_gpu, False)
                cpy_back(y.centroid_values, y.centroid_values_gpu, False)

                s.centroid_values_gpu.free()
                b.centroid_values_gpu.free()
                x.centroid_values_gpu.free()
                y.centroid_values_gpu.free()


            else:
                W1 = self.extrapolate_second_order_sw_false_block
                self.extrapolate_second_order_sw_false_func(
                    numpy.int32(N),
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
                    grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
                    )

        else:
            Domain.extrapolate_second_order_sw(self)

        if self.cotesting:
            Domain.extrapolate_second_order_sw(self.cotesting_domain)
            test_extrapolate_second_order_sw(self)
            


    # FIXME
    def ensure_numeric(A, typecode=None):
        """From numerical_tools"""
        if A is None:
            return None
        elif typecode is None:
            if isinstance(A, numpy.ndarray):
                return A
            else:
                return numpy.array(A)
        else:
            return numpy.array(A, dtype=typecode, copy=False)

    
    # FIXME
    def get_absolute(self, points):
        """From geo_reference get_absolute"""
        is_list = isinstance(poins, list)
        points = self.ensure_numeric(points, list)

        if len( points.shape) == 1:
            msg = 'Single point must have two elements'
            if not len(points) == 2:
                raise ShapeError, msg  

        msg = 'Input must be an N x 2 array or list of (x,y) values. '
        msg += 'I got an %d x %d array' %points.shape    
        if not points.shape[1] == 2:
            raise ShapeError, msg    

        if not self.mesh.geo_reference.is_absolute():
            #import copy
            #points = copy.copy(points) # Don't destroy input             
            #points[:,0] += self.mesh.geo_reference.xllcorner 
            #points[:,1] += self.mesh.geo_reference.yllcorner
            W1 = 32
            W2 = 1
            W3 = 1
            self.get_absolute_func(
                numpy.int32(points.shape[0]),
                numpy.float64(self.mesh.geo_reference.xllcorner),
                numpy.float64(self.mesh.geo_reference.yllcorner),
                self.vertex_coordinates_gpu,
                block = (W1, W2, W3),
                grid=((len(points.shape[0])+W1*W2*W3-1)/(W1*W2*W3),1)
                )


        if is_list:
            points = points.tolist()

        return points


    # FIXME
    def get_vertex_coordinates(self, triangle_id=None, absolute=False):
        if self.using_gpu:
            V = self.vertex_coordinates
            if triangle_id is None:
                if absolute is True:
                    if not self.mesh.geo_reference.is_absolute():
                        #V = self.mesh.geo_reference.get_absolute(V)
                        V = self.get_absolute(V)
                return V
            else:
                i = triangle_id
                msg =  'triangle_id must be an integer'
                assert int(i) == i, msg
                assert 0 <= i < self.number_of_triangles

                i3 = 3*i
                if absolute is True and \
                    not self.geo_reference.is_absolute():
                    
                    offset=numpy.array([self.geo_reference.get_xllcorner(),
                        self.geo_reference.get_yllcorner()], numpy.float)

                    return V[i3:i3+3,:] + offset   
                else:
                    return V[i3:i3+3,:]
        else:
            return Domain.get_vertex_coordinates(self)
                    



    # 4th level cotesting
    # Using Stream
    def manning_friction_implicit(self):
        """From shallow_water_domain"""  
        if self.using_gpu:
            N = self.number_of_elements
            W2 = 1
            W3 = 1
            
            # FIXME:
            x = self.get_vertex_coordinates()
            #cpy_back(self.vertex_coordinates, self.vertex_coordinates_gpu)
            #if (x != self.vertex_coordinates).all():
            #    print "Error: vertex_coordinates not correct"
            #    raise Exception()
                

            if self.use_sloped_mannings:
                W1 = self.manning_friction_sloped_blcok
                self.manning_friction_sloped_func.prepared_call(
                   ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                   (W1, W2, W3),
                   #self.stream[ manning_sloped_stream],
                   N,
                   self.g,
                   self.minimum_allowed_height,
   
                   self.vertex_coordinates_gpu,
                   self.quantities['stage'].centroid_values_gpu,
                   self.quantities['elevation'].vertex_values_gpu,
                   self.quantities['xmomentum'].centroid_values_gpu,
                   self.quantities['ymomentum'].centroid_values_gpu,
   
                   self.quantities['friction'].centroid_values_gpu,
                   self.quantities['xmomentum'].semi_implicit_update_gpu,
                   self.quantities['ymomentum'].semi_implicit_update_gpu
                   )
            else:
                s = self.quantities['stage']
                b = self.quantities['elevation']
                x = self.quantities['xmomentum']
                y = self.quantities['ymomentum']
                f = self.quantities['friction']
                
                s.centroid_values_gpu = get_device_array(s.centroid_values, True)
                b.vertex_values_gpu = get_device_array(b.vertex_values, True)
                x.centroid_values_gpu = get_device_array(x.centroid_values, True)
                y.centroid_values_gpu = get_device_array(y.centroid_values, True)
                f.centroid_values_gpu = get_device_array(f.centroid_values, True)
                
                x.semi_implicit_update_gpu = get_device_array(x.semi_implicit_update, True)
                y.semi_implicit_update_gpu = get_device_array(y.semi_implicit_update, True)

                W1 = self.manning_friction_flat_block
                self.manning_friction_flat_func.prepared_call(
                   ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                   (W1, W2, W3),
                   #self.stream[ manning_flat_stream ]
                   N,
                   self.g,
                   self.minimum_allowed_height,
   
                   self.quantities['stage'].centroid_values_gpu,
                   self.quantities['elevation'].vertex_values_gpu,
                   self.quantities['xmomentum'].centroid_values_gpu,
                   self.quantities['ymomentum'].centroid_values_gpu,
   
                   self.quantities['friction'].centroid_values_gpu,
                   self.quantities['xmomentum'].semi_implicit_update_gpu,
                   self.quantities['ymomentum'].semi_implicit_update_gpu
                   )
                   
                cpy_back( x.semi_implicit_update, x.semi_implicit_update_gpu, False)
                cpy_back( y.semi_implicit_update, y.semi_implicit_update_gpu, False)
                s.centroid_values_gpu.free()
                b.vertex_values_gpu.free()
                x.centroid_values_gpu.free()
                y.centroid_values_gpu.free()
                f.centroid_values_gpu.free()
                x.semi_implicit_update_gpu.free()
                y.semi_implicit_update_gpu.free()
            
        else:
            from anuga.shallow_water.shallow_water_domain import \
                    manning_friction_implicit
            manning_friction_implicit(self)

        if self.cotesting:
            from anuga.shallow_water.shallow_water_domain import \
                    manning_friction_implicit
            manning_friction_implicit(self.cotesting_domain)
            test_manning_friction_implicit(self) 




    # 4th level cotesting
    def manning_friction_explicit(self):
        """From shallow_water_domain"""  
        if self.using_gpu:
            N = self.number_of_elements
            W2 = 1
            W3 = 1
            #FIXME
            x = self.get_vertex_coordinates()
   
            if self.use_sloped_mannings:
                W1 = self.manning_friction_sloped_blcok
                self.manning_friction_sloped_func.prepared_call(
                    ((N+W1*W2*W3-1)/(W1*W2*W3),1)
                    (W1, W2, W3),
                    N,
                    self.g,
                    self.minimum_allowed_height,

                    self.vertex_coordinates_gpu,
                    self.quantities['stage'].centroid_values_gpu,
                    self.quantities['xmomentum'].centroid_values_gpu,
                    self.quantities['ymomentum'].centroid_values_gpu,
                    self.quantities['elevation'].vertex_values_gpu,

                    self.quantities['friction'].centroid_values_gpu,
                    self.quantities['xmomentum'].explicit_update_gpu,
                    self.quantities['ymomentum'].explicit_update_gpu
                    )
            else:
                W1 = self.manning_friction_flat_block
                self.manning_friction_flat_func.prepared_call(
                    ((N+W1*W2*W3-1)/(W1*W2*W3),1)
                    (W1, W2, W3),
                    N,
                    self.g,
                    self.minimum_allowed_height,

                    self.quantities['stage'].centroid_values_gpu,
                    self.quantities['xmomentum'].centroid_values_gpu,
                    self.quantities['ymomentum'].centroid_values_gpu,
                    self.quantities['elevation'].vertex_values_gpu,

                    self.quantities['friction'].centroid_values_gpu,
                    self.quantities['xmomentum'].explicit_update_gpu,
                    self.quantities['ymomentum'].explicit_update_gpu
                    )
   
        else:
            from anuga.shallow_water.shallow_water_domain import \
                    manning_friction_implicit
            manning_friction_implicit(self)



    # 3rd level cotesting
    def compute_forcing_terms(self):
        if self.using_gpu:
            for f in self.forcing_terms:
                f()
        else:
            Domain.compute_forcing_terms(self)
        if self.cotesting:
            #Domain.compute_forcing_terms(self.cotesting_domain)
            test_compute_forcing_terms(self)





    # 3rd level cotesting
    # Using Stream
    def update_conserved_quantities(self):
        if self.using_gpu :
            N = self.number_of_elements
            W1 = self.update_block 
            W2 = 1
            W3 = 1
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                Q.centroid_values_gpu = get_device_array(Q.centroid_values, True)
                Q.explicit_update_gpu = get_device_array(Q.explicit_update, True)
                Q.semi_implicit_update_gpu=get_device_array(Q.semi_implicit_update, True)
                self.update_func.prepared_call(
                    ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                    (W1, W2, W3),
                    N,
                    self.timestep,
                    Q.centroid_values_gpu,
                    Q.explicit_update_gpu,
                    Q.semi_implicit_update_gpu
                    )

                #drv.memset_d32(Q.semi_implicit_update_gpu, 0, N*2)
                Q.semi_implicit_update[:] = 0.0

                cpy_back( Q.centroid_values, Q.centroid_values_gpu, False)
                cpy_back( Q.semi_implicit_update,Q.semi_implicit_update_gpu,False)
                Q.centroid_values_gpu.free()
                Q.explicit_update_gpu.free()
                Q.semi_implicit_update_gpu.free()
        else:
            Domain.update_conserved_quantities(self)

        if self.cotesting:
            Domain.update_conserved_quantities(self.cotesting_domain)
            test_update_conserved_quantities(self)




    # 3rd level cotesting
    # Using asynchronous_transfer
    def backup_conserved_quantities(self):
        if self.using_gpu:
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                #FIXME: may use asynchronous_transfer
                #drv.memcpy_dtod(
                #        Q.centroid_backup_values_gpu, 
                #        Q.centroid_values_gpu)
                drv.memcpy_dtod_async(
                        Q.centroid_backup_values_gpu, 
                        Q.centroid_values_gpu,
                        size = Q.centroid_values.nbytes,
                        stream = drv.Stream()
                        )
                        
        else:
            Domain.backup_conserved_quantities(self)

        if self.cotesting:
            Domain.backup_conserved_quantities(self.cotesting_domain)




    # 3rd level cotesting
    # Using Stream
    def saxpy_conserved_quantities(self, a, b):
        if self.using_gpu:
            N = self.number_of_elements
            W1 = self.saxpy_centroid_values_block
            W2 = 1
            W3 = 1
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                self.saxpy_centroid_values_func.prepared_call(
                    ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                    (W1, W2, W3),
                    # FIXME: use 3 different stream 
                    #self.stream[ saxpy_stream ],
                    N,
                    a,
                    b,
                    Q.centroid_values_gpu,
                    Q.centroid_backup_values_gpu
                    )
                    
        else:
            Domain.saxpy_conserved_quantities(self, a, b)


        if self.cotesting:
            Domain.saxpy_conserved_quantities(self.cotesting_domain, a, b)
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                c1 = Q.centroid_values
                cpy_back( c1, Q.centroid_values_gpu)
                c2= self.cotesting_domain.quantities[name].centroid_values
                if not numpy.allclose(c1, c2):
                    print "Error: saxpy_conserved_quantities", name

                
                

                
    # 3rd level 
    # Using Stream
    def update_centroids_of_velocities_and_height(self):
        if self.using_gpu:
            N = self.number_of_elements
            Nb = self.quantities['stage'].boundary_values.shape[0]
            W1 = self.set_boundary_block
            W2 = 1
            W3 = 1

            s = self.quantities['stage']
            b = self.quantities['elevation']
            xm = self.quantities['xmomentum']
            ym = self.quantities['ymomentum']
            h = self.quantities['height']
            xv = self.quantities['xvelocity']
            yv = self.quantities['yvelocity']

            self.boundary_cells_gpu = get_device_array(self.boundary_cells, True)
            self.boundary_edges_gpu = get_device_array(self.boundary_edges, True)
            b.boundary_values_gpu = get_device_array(b.boundary_values, True)
            b.edge_values_gpu = get_device_array(b.edge_values, True)

            self.set_boundary_values_from_edges_func.prepared_call(
                ((Nb+W1*W2*W3-1)/(W1*W2*W3),1),
                (W1, W2, W3),
                #self.stream[ set_boundary_from_edge_stream],
                Nb,
                self.boundary_cells_gpu,
                self.boundary_edges_gpu,
                self.quantities['elevation'].boundary_values_gpu,
                self.quantities['elevation'].edge_values_gpu
                )

            cpy_back( b.boundary_values, b.boundary_values_gpu, False)
            self.boundary_cells_gpu.free()
            self.boundary_edges_gpu.free()
            b.boundary_values_gpu.free()
            b.edge_values_gpu.free()

            s.centroid_values_gpu = get_device_array(s.centroid_values, True)
            xm.centroid_values_gpu = get_device_array(xm.centroid_values, True)
            ym.centroid_values_gpu = get_device_array(ym.centroid_values, True)
            h.centroid_values_gpu = get_device_array(h.centroid_values, True)
            b.centroid_values_gpu = get_device_array(b.centroid_values, True)
            xv.centroid_values_gpu = get_device_array(xv.centroid_values, True)
            yv.centroid_values_gpu = get_device_array(yv.centroid_values, True)

            s.boundary_values_gpu = get_device_array(s.boundary_values, True)
            xm.boundary_values_gpu = get_device_array(xm.boundary_values, True)
            ym.boundary_values_gpu = get_device_array(ym.boundary_values, True)
            h.boundary_values_gpu = get_device_array(h.boundary_values, True)
            b.boundary_values_gpu = get_device_array(b.boundary_values, True)
            xv.boundary_values_gpu = get_device_array(xv.boundary_values, True)
            yv.boundary_values_gpu = get_device_array(yv.boundary_values, True)

            W1 = self.update_centroids_of_velocities_and_height_block
            self.update_centroids_of_velocities_and_height_func.\
                prepared_call(
                    ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                    (W1, W2, W3),
                    #self.stream[ set_boundary_from_edge_stream],
                    N,
                    Nb,
                    self.quantities['stage'].centroid_values_gpu,
                    self.quantities['xmomentum'].centroid_values_gpu,
                    self.quantities['ymomentum'].centroid_values_gpu,
                    self.quantities['height'].centroid_values_gpu,
                    self.quantities['elevation'].centroid_values_gpu,
                    self.quantities['xvelocity'].centroid_values_gpu,
                    self.quantities['yvelocity'].centroid_values_gpu,

                    self.quantities['stage'].boundary_values_gpu,
                    self.quantities['xmomentum'].boundary_values_gpu,
                    self.quantities['ymomentum'].boundary_values_gpu,
                    self.quantities['height'].boundary_values_gpu,
                    self.quantities['elevation'].boundary_values_gpu,
                    self.quantities['xvelocity'].boundary_values_gpu,
                    self.quantities['yvelocity'].boundary_values_gpu
                )
            
            cpy_back( h.centroid_values, h.centroid_values_gpu, False)
            cpy_back( xv.centroid_values,xv.centroid_values_gpu, False)
            cpy_back( yv.centroid_values,yv.centroid_values_gpu, False)
            cpy_back( h.boundary_values, h.boundary_values_gpu, False)
            cpy_back( xv.boundary_values, xv.boundary_values_gpu, False)
            cpy_back( yv.boundary_values, yv.boundary_values_gpu, False)

            s.centroid_values_gpu.free() 
            xm.centroid_values_gpu.free()
            ym.centroid_values_gpu.free()
            h.centroid_values_gpu.free() 
            b.centroid_values_gpu.free() 
            xv.centroid_values_gpu.free()
            yv.centroid_values_gpu.free()

            s.boundary_values_gpu.free() 
            xm.boundary_values_gpu.free()
            ym.boundary_values_gpu.free()
            h.boundary_values_gpu.free() 
            b.boundary_values_gpu.free() 
            xv.boundary_values_gpu.free()
            yv.boundary_values_gpu.free()

        else:
            Domain.update_centroids_of_velocities_and_height(self)


        if self.cotesting:
            Domain.update_centroids_of_velocities_and_height(
                    self.cotesting_domain)
            test_update_centroids_of_velocities_and_height(self)



    

    # 3th level cotesting
    def compute_fluxes(self):
        if self.using_gpu :
            N = self.number_of_elements
            W2 = 1
            W3 = 1
            

            if self.compute_fluxes_method == 'original':
                print "original"
            elif self.compute_fluxes_method == 'wb_1':
                print "wb_1"
            elif self.compute_fluxes_method == 'wb_2':
                s = self.quantities['stage']
                b = self.quantities['elevation']
                x = self.quantities['xmomentum']
                y = self.quantities['ymomentum']

                self.timestep_array_gpu = get_device_array(self.timestep_array, True)
                self.neighbours_gpu = get_device_array(self.neighbours, True)
                self.neighbour_edges_gpu = get_device_array(self.neighbour_edges, True)
                self.normals_gpu = get_device_array(self.normals, True)
                self.edgelengths_gpu = get_device_array(self.edgelengths, True)
                self.radii_gpu = get_device_array(self.radii, True)
                self.areas_gpu = get_device_array(self.areas, True)
                self.tri_full_flag_gpu = get_device_array(self.tri_full_flag, True)
                self.max_speed_gpu = get_device_array(self.max_speed, True)
                self.vertex_coordinates_gpu = get_device_array(self.vertex_coordinates, True)

                s.vertex_values_gpu=get_device_array(s.vertex_values, True)
                s.centroid_values_gpu=get_device_array(s.centroid_values, True)
                b.centroid_values_gpu=get_device_array(b.centroid_values, True)

                s.edge_values_gpu=get_device_array(s.edge_values, True)
                b.edge_values_gpu=get_device_array(b.edge_values, True)
                x.edge_values_gpu=get_device_array(x.edge_values, True)
                y.edge_values_gpu=get_device_array(y.edge_values, True)


                s.boundary_values_gpu=get_device_array(s.boundary_values, True)
                x.boundary_values_gpu=get_device_array(x.boundary_values, True)
                y.boundary_values_gpu=get_device_array(y.boundary_values, True)
            
                s.explicit_update_gpu=get_device_array(s.explicit_update, True)
                x.explicit_update_gpu=get_device_array(x.explicit_update, True)
                y.explicit_update_gpu=get_device_array(y.explicit_update, True)
            
                W1 = self.compute_fluxes_central_structure_block
                self.compute_fluxes_func.prepared_call(
                        ((N+W1*W2*W3-1)/(W1*W2*W3), 1),
                        (W1, W2, W3),
                        #self.stream[ cf_central_stream ],
                        N,
                        self.evolve_max_timestep,
                        self.g,
                        self.epsilon,
                        self.H0 * self.H0,
                        self.H0 * 10,
                        self.optimise_dry_cells,
                        
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
                        self.max_speed_gpu
                        )
                    
                #drv.memcpy_dtoh( self.timestep_array, 
                #        self.timestep_array_gpu)
                drv.memcpy_dtoh_async( self.timestep_array, 
                        self.timestep_array_gpu, 
                        stream = self.stream[ cf_central_stream]
                        ) 

                W1 = self.gravity_wb_block
                self.gravity_wb_func.prepared_call(
                        ((N + W1*W2*W3-1)/(W1*W2*W3),1),
                        (W1, W2, W3),
                        #self.stream[ cf_central_stream ],
                        N,
                        self.g,
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
                        self.edgelengths_gpu
                        )

                cpy_back( self.max_speed, self.max_speed_gpu, False)
                cpy_back( s.explicit_update, s.explicit_update_gpu, False)
                cpy_back( x.explicit_update, x.explicit_update_gpu, False)
                cpy_back( y.explicit_update, y.explicit_update_gpu, False)

                self.timestep_array_gpu.free()
                self.neighbours_gpu.free()
                self.neighbour_edges_gpu.free()
                self.normals_gpu.free()
                self.edgelengths_gpu.free()
                self.radii_gpu.free()
                self.areas_gpu.free()
                self.tri_full_flag_gpu.free()
                self.max_speed_gpu.free()
                self.vertex_coordinates_gpu.free()

                s.edge_values_gpu.free()
                x.edge_values_gpu.free()
                y.edge_values_gpu.free()
                b.edge_values_gpu.free()

                s.boundary_values_gpu.free()
                x.boundary_values_gpu.free()
                y.boundary_values_gpu.free()

                s.explicit_update_gpu.free()
                x.explicit_update_gpu.free()
                y.explicit_update_gpu.free()

            elif self.compute_fluxes_method == 'wb_3':
                print "wb_3"
            elif self.compute_fluxes_method == 'tsunami':
                print "tsunami"
            else:
                raise Exception('unknown compute_fluxes_method')

            
            
            self.flux_timestep = numpy.min(self.timestep_array)
            #print self.flux_timestep
            #print numpy.min(self.timestep_array)
            

        else:
            Domain.compute_fluxes(self)
            #print self.flux_timestep


        if self.cotesting:
            Domain.compute_fluxes(self.cotesting_domain)
            test_compute_fluxes(self)




    # For cotesting purpose
    # 3rd level cotesting purpose
    def update_timestep(self, yieldstep, finaltime):
        Domain.update_timestep(self, yieldstep, finaltime)
        if self.cotesting:
            Domain.update_timestep(
                    self.cotesting_domain, yieldstep, finaltime)
            
            test_update_timestep(self)
            
        


    # For cotesting purpose
    # 2nd level cotesting
    def apply_fractional_steps(self):
        pass
        #Domain.apply_fractional_steps(self)
        #if self.cotesting:
        #    Domain.apply_fractional_steps(self.cotesting_domain)
            


    # For cotesting purpose
    # 2nd level cotesting
    def update_ghosts(self):
        Domain.update_ghosts(self)
        if self.cotesting:
            Domain.update_ghosts(self.cotesting_domain)
            test_update_ghosts(self)



    # 2nd level cotesting
    def distribute_to_vertices_and_edges(self):
        if  self.using_gpu:
            N = self.number_of_elements
            W2 = 1
            W3 = 1

            if self.compute_fluxes_method == 'tsunami':
                W1 = self.protect_swb2_block
                self.protect_swb2_func.prepared_call(
                        ((N+W1*W2*W3-1)/(W1*W2*W3),1)
                        (W1, W2, W3),
                        N,
                        self.minimum_allowed_height,
                        self.maximum_allowed_speed,
                        number_of_elements.float64(self.epsilon),
                        self.quantities['stage'].centroid_values_gpu,
                        self.quantities['stage'].vertex_values_gpu,
                        self.quantities['elevation'].centroid_values_gpu,
                        self.quantities['elevation'].vertex_values_gpu,
                        self.quantities['xmomentum'].centroid_values_gpu,
                        self.quantities['ymomentum'].centroid_values_gpu,
                        self.areas_gpu
                        )

                W1 = self.extrapolate_second_order_edge_swb2_block
                self.extrapolate_second_order_edge_swb2_func(
                        ((N+W1*W2*W3-1)/(W1*W2*W3),1)
                        (W1, W2, W3),
                        N,
                        self.optimise_dry_cells,
                        self.extrapolate_velocity_second_order,
                        self.minimum_allowed_height,
                        self.beta_w,
                        self.beta_w_dry,
                        self.beta_uh,
                        self.beta_uh_dry,
                        self.beta_vh,
                        self.beta_vh_dry,

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
                        self.mim_elevation_edgevalue_gpu,
                        self.max_elevation_edgevalue_gpu,
                        self.count_wet_neighbours_gpu
                        )
                
            elif self.use_edge_limiter:
                self.protect_against_infinitesimal_and_negative_heights()

                for name in self.conserved_quantities:
                    Q = self.quantities[name]
                    if self._order_ == 1:
                        W1 = self.extrapolate_first_order_block
                        #Q.extrapolate_first_order()
                        self.extrapolate_first_order_func.prepared_call(
                            ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                            (W1, W2, W3),
                            N,
                            Q.centroid_values_gpu,
                            Q.edge_values_gpu,
                            Q.vertex_values_gpu
                            )

                        drv.memset_d32(Q.x_gradient_gpu,0,N*2)
                        drv.memset_d32(Q.y_gradient_gpu,0,N*2)
                        

                    elif self._order_ == 2:
                        #Q.extrapolate_second_order_and_limit_by_edge()
                        W1 = self.\
                        extrapolate_second_order_and_limit_by_edge_block
                        self.\
                        extrapolate_second_order_and_limit_by_edge_func.\
                        prepared_call(
                                ((N+W1*W2*W3-1)/(W1*W2),1),
                                (W1, W2, W3),
                                N,
                                Q.beta,
                                self.centroid_coordinates_gpu,
                                self.vertex_coordinates_gpu,
                                self.number_of_boundaries_gpu,
                                self.surrogate_neighbours_gpu,
                                self.neighbours_gpu,

                                Q.centroid_values_gpu,
                                Q.vertex_values_gpu,
                                Q.edge_values_gpu,
                                Q.x_gradient_gpu,
                                Q.y_gradient_gpu
                                )
                    else:
                        raise Exception('Unknown order')


                self.balance_deep_and_shallow()

                for name in self.conserved_quantities:
                    Q = self.quantities[name]
                    N = Q.vertex_values.shape[0]
                    #Q.interpolate_from_vertices_to_edges()
                    W1 = self.interpolate_from_vertices_to_edges_block
                    self.interpolate_from_vertices_to_edges_func.\
                            prepared_call(
                                    ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                                    (W1, W2, W3),
                                    N,
                                    Q.vertex_values_gpu,
                                    Q.edge_values_gpu
                                    )

            # using vertex limiter
            else:
                self.protect_against_infinitesimal_and_negative_heights()

                if self.optimised_gradient_limiter:
                    if self._order_ == 1:
                        for name in self.conserved_quantities:
                            Q = self.quantities[name]
                            W1 = self.extrapolate_first_order_block
                            self.extrapolate_first_order_func.\
                                    prepared_call(
                                        ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                                        (W1, W2, W3),
                                        N,
                                        Q.centroid_values_gpu,
                                        Q.edge_values_gpu,
                                        Q.vertex_values_gpu,
                                        )

                            drv.memset_d32( Q.x_gradient_gpu, 0, N*2)
                            drv.memset_d32( Q.y_gradient_gpu, 0, N*2)
                            
                            # cotesting point
                            if self.cotesting:
                                Q2 = self.cotesting_domain.quantities[name]
                                Q.extrapolate_first_order()
                                test_extrapolate_first_order(self)

                    elif self._order_ == 2:
                        self.extrapolate_second_order_sw()
                    else:
                        raise Exception('Unknown order')

                else:
                    for name in self.conserved_quantities:
                        Q = self.quantities[name]

                        if self._order_ == 1:
                            #Q.extrapolate_first_order()
                            W1 = self.extrapolate_first_order_block
                            self.extrapolate_first_order_func.\
                                    prepared_call(
                                            ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                                            (W1, W2, W3),
                                            N,
                                            Q.centroid_values_gpu,
                                            Q.edge_values_gpu,
                                            Q.vertex_values_gpu,
                                            )

                            drv.memset_d32( Q.x_gradient_gpu, 0, N*2)
                            drv.memset_d32( Q.y_gradient_gpu, 0, N*2)


                            # cotesting point
                            if self.cotesting:
                                Q2 = self.cotesting_domain.quantities[name]
                                Q2.extrapolate_first_order()
                                test_extrapolate_first_order(self)

                        elif self._order_ == 2:
                            #Q.extrapolate_second_order_and_limit_by_vertex
                            W1 = self.extrapolate_second_order_and_limit_by_vertex_block
                            self.extrapolate_second_order_and_limit_by_vertex_func.prepared_call(
                                ((N+W1*W2*W3-1)/(W1*W2),1),
                                (W1, W2, W3),
                                N,
                                Q.beta,
                                self.centroid_coordinates_gpu,
                                self.vertex_coordinates_gpu,
                                self.number_of_boundaries_gpu,
                                self.surrogate_neighbours_gpu,
                                self.neighbours_gpu,

                                Q.centroid_values_gpu,
                                Q.vertex_values_gpu,
                                Q.edge_values_gpu,
                                Q.x_gradient_gpu,
                                Q.y_gradient_gpu
                                )
                                
                            # cotesting point
                            if self.cotesting:
                                Q2 = self.cotesting_domain.quantities[name]
                                Q2.extrapolate_second_order_and_limit_by_vertex()
                                test_extrapolate_second_order_and_limit_by_vertex(self)

                        else:
                            raise Exception('Unknown order')

                self.balance_deep_and_shallow()

                # cotesting point input of interpolate_from_vertices_to_edges
                if self.cotesting:
                    test_interpolate_from_vertices_to_edges(self)
                    
                for name in self.conserved_quantities:
                    Q = self.quantities[name]
                    
                    Q.vertex_values_gpu = get_device_array(Q.vertex_values, True) 
                    Q.edge_values_gpu = get_device_array(Q.edge_values, True) 

                    W1 = self.interpolate_from_vertices_to_edges_block
                    self.interpolate_from_vertices_to_edges_func.\
                            prepared_call(
                                    ( (N+W1*W2*W3-1)/(W1*W2*W3) ,1),
                                    (W1, W2, W3),
                                    N,
                                    Q.vertex_values_gpu,
                                    Q.edge_values_gpu,
                                    )
                    cpy_back( Q.vertex_values, Q.vertex_values_gpu, False)
                    cpy_back( Q.edge_values, Q.edge_values_gpu, False)
                    
                    Q.vertex_values_gpu.free()
                    Q.edge_values_gpu.free()

                ## cotesting point output of interpolate_from_vertices_to_edges
                    if self.cotesting:
                        Q2 = self.cotesting_domain.quantities[name]
                        Q2.interpolate_from_vertices_to_edges()
                        test_interpolate_from_vertices_to_edges(self)

        else:
            Domain.distribute_to_vertices_and_edges(self)
            if self.cotesting:
                self.cotesting_domain.distribute_to_vertices_and_edges()


        if self.cotesting:
            test_distribute_to_vertexs_and_edges(self)



    # 2nd level cotesting
    def update_boundary(self):
        if self.using_gpu:
            W2 = 1
            W3 = 1

            for tag in self.tag_boundary_cells:
                B = self.boundary_map[tag]
                if B is None:
                    continue
                
                #segment_edges = self.tag_boundary_cells[tag]
                ids = self.tag_boundary_cells[tag]
                N = len(ids)

                # def evaluate_segment(self, domain, segment_edges):
                if ids is None:
                    continue
                s = self.quantities['stage']
                b = self.quantities['elevation']
                h = self.quantities['height']
                xm = self.quantities['xmomentum']
                ym = self.quantities['ymomentum']
                xv = self.quantities['xvelocity']
                yv = self.quantities['yvelocity']

                self.boundary_index = []
                self.boundary_index.append( get_device_array(numpy.asarray(ids), True))
                self.boundary_index.append( get_device_array( numpy.asarray(self.boundary_cells[ids]), True))
                self.boundary_index.append( get_device_array( numpy.asarray(self.boundary_edges[ids]), True))
                self.normals_gpu = get_device_array(self.normals, True)
                
                s.edge_values_gpu = get_device_array(s.edge_values, True)
                b.edge_values_gpu = get_device_array(b.edge_values, True)
                xm.edge_values_gpu = get_device_array(xm.edge_values, True)
                ym.edge_values_gpu = get_device_array(ym.edge_values, True)
                h.edge_values_gpu = get_device_array(h.edge_values, True)
                xv.edge_values_gpu = get_device_array(xv.edge_values, True)
                yv.edge_values_gpu = get_device_array(yv.edge_values, True)

                s.boundary_values_gpu = get_device_array(s.boundary_values, True)
                b.boundary_values_gpu = get_device_array(b.boundary_values, True)
                xm.boundary_values_gpu =get_device_array(xm.boundary_values, True)
                ym.boundary_values_gpu =get_device_array(ym.boundary_values, True)
                h.boundary_values_gpu = get_device_array(h.boundary_values, True)
                xv.boundary_values_gpu =get_device_array(xv.boundary_values, True)
                yv.boundary_values_gpu =get_device_array(yv.boundary_values, True)

                if isinstance(B, Reflective_boundary):
                    W1 = self.evaluate_segment_reflective_block
                    self.evaluate_segment_reflective_func.prepared_call(
                        ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                        (W1, W2, W3),
                        self.number_of_elements,
                        N,
                        self.boundary_index[0],
                        self.boundary_index[1],
                        self.boundary_index[2],

                        self.normals_gpu,
                        self.quantities['stage'].edge_values_gpu,
                        self.quantities['elevation'].edge_values_gpu,
                        self.quantities['height'].edge_values_gpu,
                        self.quantities['xmomentum'].edge_values_gpu,
                        self.quantities['ymomentum'].edge_values_gpu,
                        self.quantities['xvelocity'].edge_values_gpu,
                        self.quantities['yvelocity'].edge_values_gpu,

                        self.quantities['stage'].boundary_values_gpu,
                        self.quantities['elevation'].boundary_values_gpu,
                        self.quantities['height'].boundary_values_gpu,
                        self.quantities['xmomentum'].boundary_values_gpu,
                        self.quantities['ymomentum'].boundary_values_gpu,
                        self.quantities['xvelocity'].boundary_values_gpu,
                        self.quantities['yvelocity'].boundary_values_gpu
                        )
                    
                    cpy_back( s.boundary_values, s.boundary_values_gpu, False)   
                    cpy_back( b.boundary_values, b.boundary_values_gpu, False)   
                    cpy_back( xm.boundary_values, xm.boundary_values_gpu, False)   
                    cpy_back( ym.boundary_values, ym.boundary_values_gpu, False)   
                    cpy_back( h.boundary_values, h.boundary_values_gpu, False)   
                    cpy_back( xv.boundary_values, xv.boundary_values_gpu, False)   
                    cpy_back( yv.boundary_values, yv.boundary_values_gpu, False)   

                    s.boundary_values_gpu.free() 
                    b.boundary_values_gpu.free()
                    xm.boundary_values_gpu.free()
                    ym.boundary_values_gpu.free()
                    h.boundary_values_gpu.free() 
                    xv.boundary_values_gpu.free()
                    yv.boundary_values_gpu.free()

                elif isinstance(B, Dirichlet_boundary):
                    q_bdry = B.dirichlet_values
                    conserved_quantities = True
                    if len(q_bdry) == len(self.evolved_quantities):
                        conserved_quantities = False
                    if  conserved_quantities:
                        for j, name in enumerate(self.evolved_quantities):
                            Q = self.quantities[name]
                            W1 = self.evaluate_segment_dirichlet_1_block
                            self.evaluate_segment_dirichlet_1_func.\
                                prepared_call(
                                        ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                                        (W1, W2, W3),
                                        N,
                                        self.boundary_index[tag][0],
                                        self.boundary_index[tag][1],
                                        self.boundary_index[tag][2],

                                        Q.boundary_values_gpu,
                                        Q.edge_values_gpu
                                        )


                    if conserved_quantities:
                        quantities = self.conserved_quantities
                    else:
                        quantities = self.evolved_quantities

                    for j, name in enumerate(quantities):
                        Q = self.quantities[name]
                        W1 = self.evaluate_segment_dirichlet_2_block
                        self.evaluate_segment_dirichlet_2_func.\
                                prepared_call(
                                        ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                                        (W1, W2, W3),
                                        N,
                                        q_bdry[j],
                                        self.boundary_index[tag][0],
                                        Q.boundary_values_gpu
                                        )

                else:
                    raise Exception("Can not find right type")
            
        else:
            Generic_Domain.update_boundary(self)


        if self.cotesting:
            Generic_Domain.update_boundary(self.cotesting_domain)
            test_update_boundary(self)



    # 2nd level cotesting
    # Using Stream
    def update_other_quantities(self):
        if self.using_gpu:
            if self.flow_algorithm == 'yusuke':
                return
            
            
            self.update_centroids_of_velocities_and_height()
            
            N = self.number_of_elements
            W1 = self.extrapolate_first_order_block
            W2 = 1
            W3 = 1
            for name in ['height', 'xvelocity', 'yvelocity']:
                Q = self.quantities[name]
                Q.centroid_values_gpu = get_device_array(Q.centroid_values, True)
                Q.edge_values_gpu = get_device_array(Q.edge_values, True)
                Q.vertex_values_gpu=get_device_array(Q.vertex_values, True)

                self.extrapolate_first_order_func.prepared_call(
                        ((N+W1*W2*W3-1)/(W1*W2*W3),1),
                        (W1, W2, W3),
                        N,
                        Q.centroid_values_gpu,
                        Q.edge_values_gpu,
                        Q.vertex_values_gpu
                        )

                #drv.memset_d32(Q.x_gradient_gpu,0,N*2)
                #drv.memset_d32(Q.y_gradient_gpu,0,N*2)
                Q.x_gradient[:] = 0
                Q.y_gradient[:] = 0
                cpy_back(Q.edge_values, Q.edge_values_gpu, False)
                cpy_back(Q.vertex_values, Q.vertex_values_gpu, False)
                Q.centroid_values_gpu.free()
                Q.edge_values_gpu.free()
                Q.vertex_values_gpu.free()
        else:
            Domain.update_other_quantities(self)


        if self.cotesting:
            for name in ['height', 'xvelocity', 'yvelocity']:
                Q = self.cotesting_domain.quantities[name]
                Q.extrapolate_first_order()
            test_update_other_quantities(self)



    # For cotesting purpose
    # 2nd level cotesting
    def update_extrema(self):
        Domain.update_extrema(self)
        if self.cotesting:
            Domain.update_extrema(self.cotesting_domain)
            test_update_extrema(self)



    # For cotesting purpose
    # 2nd level cotesting
    def store_timestep(self):
        if self.using_gpu:
            self.copy_back_necessary_data()
        self.writer.store_timestep()
        if self.cotesting_domain:
            self.cotesting_domain.writer.store_timestep()

    
    
    # For cotesting purpose
    # 2nd level cotesting
    def evolve_one_euler_step(self, yieldstep, finaltime):
        Domain.evolve_one_euler_step(self, yieldstep, finaltime)
        if self.cotesting:
            #Domain.evolve_one_euler_step(self.cotesting_domain, 
            #        yieldstep, finaltime)
            test_evolve_one_euler_step(self)




    # For cotesting purpose
    # 2nd level cotesting
    def evolve_one_rk2_step(self, yieldstep, finaltime):
        Domain.evolve_one_rk2_step(self, yieldstep, finaltime)
        if self.cotesting:
            Domain.evolve_one_rk2_step(self.cotesting_domain, 
                    yieldstep, finaltime)
            s1 = self.quantities['stage']
            x1 = self.quantities['xmomentum']
            y1 = self.quantities['ymomentum']
            e1 = self.quantities['elevation']

            cpy_back( s1.centroid_values, s1.centroid_values_gpu)
            cpy_back( s1.vertex_values, s1.vertex_values_gpu)
            cpy_back( x1.centroid_values, x1.centroid_values_gpu)
            cpy_back( y1.centroid_values, y1.centroid_values_gpu)
            cpy_back( e1.centroid_values, e1.centroid_values_gpu)
            cpy_back( e1.vertex_values, e1.vertex_values_gpu)


            s2 = self.cotesting_domain.quantities['stage']
            x2 = self.cotesting_domain.quantities['xmomentum']
            y2 = self.cotesting_domain.quantities['ymomentum']
            e2 = self.cotesting_domain.quantities['elevation']

                
            ipt = []
            ipt.append( numpy.allclose( s1.centroid_values, 
                        s2.centroid_values))
            ipt.append( numpy.allclose( s1.vertex_values, 
                        s2.vertex_values))
            ipt.append( numpy.allclose( x1.centroid_values, 
                        x2.centroid_values))
            ipt.append( numpy.allclose( y1.centroid_values, 
                        y2.centroid_values))
            ipt.append( numpy.allclose( e1.centroid_values,
                        e2.centroid_values))
            ipt.append( numpy.allclose( e1.vertex_values,
                        e2.vertex_values))

            if not ipt.count(True) == ipt.__len__():
                print "   --> evolve_one_rk2_step ", ipt



    # For cotesting purpose
    # 2nd level cotesting
    def evolve_one_rk3_step(self, yieldstep, finaltime):
        Domain.evolve_one_rk3_step(self, yieldstep, finaltime)
        if self.cotesting:
            Domain.evolve_one_rk3_step(self.cotesting_domain, 
                    yieldstep, finaltime)
            s1 = self.quantities['stage']
            x1 = self.quantities['xmomentum']
            y1 = self.quantities['ymomentum']
            e1 = self.quantities['elevation']

            cpy_back( s1.centroid_values, s1.centroid_values_gpu)
            cpy_back( s1.vertex_values, s1.vertex_values_gpu)
            cpy_back( x1.centroid_values, x1.centroid_values_gpu)
            cpy_back( y1.centroid_values, y1.centroid_values_gpu)
            cpy_back( e1.centroid_values, e1.centroid_values_gpu)
            cpy_back( e1.vertex_values, e1.vertex_values_gpu)


            s2 = self.cotesting_domain.quantities['stage']
            x2 = self.cotesting_domain.quantities['xmomentum']
            y2 = self.cotesting_domain.quantities['ymomentum']
            e2 = self.cotesting_domain.quantities['elevation']

                
            ipt = []
            ipt.append( numpy.allclose( s1.centroid_values, 
                        s2.centroid_values))
            ipt.append( numpy.allclose( s1.vertex_values, 
                        s2.vertex_values))
            ipt.append( numpy.allclose( x1.centroid_values, 
                        x2.centroid_values))
            ipt.append( numpy.allclose( y1.centroid_values, 
                        y2.centroid_values))
            ipt.append( numpy.allclose( e1.centroid_values,
                        e2.centroid_values))
            ipt.append( numpy.allclose( e1.vertex_values,
                        e2.vertex_values))

            if not ipt.count(True) == ipt.__len__():
                print "   --> evolve_one_rk3_step ", ipt


    # 1st level cotesting
    def evolve(self, 
                yieldstep=None,
                finaltime=None,
                duration=None,
                skip_initial_step=False):
        print " --> Number of elements: %d" % self.number_of_elements


        if self.using_gpu:
            """ Prepare to use GPU version calculating evolve procedure
            """

            print " *** Enable GPU Evolve ***"
            from anuga_cuda import sort_domain
            sort_domain(self)
            

            
            """ This is for testing purpose.
                
            We shadow copy current generated domain variable, and 
            proceed same evolve procedure with ANUGA original functions
            step by step, and compare the result generated by each step 
            to check whether they agree on each other.
            """

            if self.cotesting:
                print " *** Enable Cotesting ***"


                import copy
                self.cotesting_domain = copy.deepcopy(self)
                self.cotesting_domain.using_gpu = False
                self.cotesting_domain.cotesting = False



            """ Rearranged domain method is enabled here

            Basically, all the multi-dimension arrays are rearranged:
                Original:   
                            AAA
                            BBB
                            CCC
                            ...

                After rearranged:
                            ABC...
                            ABC...
                            ABC...

            so that, in most cases, threads can address the data stored nearby.
            This helps meet the requirements for coalsced memory access.
            """
            
            if self.rearranged_domain:
                print " *** Enable Rearraned Domain ***"
                from anuga_cuda import rearrange_domain
                self = rearrange_domain(self, False)
                # Since rearranged domain is defaulted to be sorted
                #if self.cotesting:
                #    sort_domain(self.cotesting_domain)



            """ CUDA stream technique is enabled here

            If device has the capabiliby and this technique is enabled,
            we need to generate a certain number of streams first. Otherwise,
            this list is filled up with None, which stands for the default 
            main stream
            """
            
            self.stream = []
            if self.using_stream:
                print " *** Enable Strem ***"
                for i in range(kbc.__len__()):
                    self.stream.append(drv.Stream())
            else:
                for i in range(kbc.__len__()):
                    self.stream.append(None)
            
            

            """ Fix forcing_terms 

            'forcing_terms' is a list of function pointers. 
            """
            
            from anuga.shallow_water.shallow_water_domain import \
                    manning_friction_implicit, manning_friction_explicit

            f = self.forcing_terms
            for i in range(len(f)):
                if f[i] == manning_friction_implicit:
                    f[i] = self.manning_friction_implicit
                elif f[i] == manning_friction_explicit:
                    f[i] == self.manning_friction_explicit
                else:
                    print "Error: in fixing forcing_terms", f[i]
                    raise Exception()
                #FIXME
                # in shallow_water_domain set_gravity_methon function
                # gravity functions can be the elements of forcing_terms

            

            """ Get kernels and transfer data to device
            """
            
            self.equip_kernel_functions()
                    
            #self.lock_array_page()
            
            
            # Initial start time
            ini_time = time.time()
            start = drv.Event()
            end = drv.Event()
            start.record()

            #self.allocate_device_array()

            #self.timestep_array[:] = self.evolve_max_timestep
            
            #self.asynchronous_transfer()
            #ctx.synchronize()
            
            # Finish data allocating and copying in
            ini_evo = time.time()
            end.record()
            end.synchronize()
            secs = start.time_till(end)*1e-3
            print "Data copy in time: %3.7f" % secs
        else:
            ini_time = time.time()
            ini_evo = time.time()
            

        
        if not self.using_gpu and self.cotesting:
            """ This is for testing purpose.
                
            We shadow copy current generated domain variable, and 
            proceed same evolve procedure with ANUGA original functions
            step by step, and compare the result generated by each step 
            to check whether they agree on each other.
            """
            print " *** Enable Cotesting ***"


            import copy
            self.cotesting_domain = copy.deepcopy(self)
            self.cotesting_domain.using_gpu = False
            self.cotesting_domain.cotesting = False



        """ Start evolve procedures

        We use the 'evolve' function from Shallow_Water_domain, to proceed 
        the evolve procedure, but necessary functions are reloaded as 
        kernel functions.
        """
            
        st_time = time.time()

        #for t in Domain.evolve(self, yieldstep=yieldstep,
        #            finaltime=finaltime,
        #            duration=duration,
        #            skip_initial_step=skip_initial_step):
        #    yield(t)


        from anuga.config import  epsilon


        # Cotesting point
        self.distribute_to_vertices_and_edges()



        if self.get_time() != self.get_starttime():
            self.set_time( self.get_starttime() )

        if yieldstep is None:
            yieldstep = self.evolve_max_timestep
        else:
            yieldstep = float(yieldstep)

        self._order_ = self.default_order

        assert finaltime >= self.get_starttime()

        if finaltime is not None and duration is not None:
            raise Exception()
        else:
            if finaltime is not None:
                self.finaltime = float( finaltime )
            if duration is not None:
                self.finaltime = self.starttime + float(duration)

        N = len(self)
        self.yieldtime = self.get_time() + yieldstep

        self.recorded_min_timestep = self.evolve_max_timestep
        self.recorded_max_timestep = self.evolve_min_timestep
        self.number_of_steps = 0
        self.number_of_first_order_steps = 0


        
        # Cotesting point
        if self.cotesting:
            sc = self.cotesting_domain

            if sc.get_time() != sc.get_starttime():
                sc.set_time( sc.get_starttime())

            sc._order_ = sc.default_order

            assert finaltime >= sc.get_starttime()

            if finaltime is not None and duration is not None:
                raise Exception()
            else:
                if finaltime is not None:
                    sc.finaltime = float(finaltime)
                if duration is not None:
                    sc.finaltime = sc.starttime + float(duration)
            sc.yieldtime = sc.get_time() + yieldstep

            sc.recorded_min_timestep = sc.evolve_max_timestep
            sc.recorded_max_timestep = sc.evolve_min_timestep
            sc.number_of_steps = 0
            sc.number_of_first_order_steps = 0





        # Cotesting point
        self.update_ghosts()
        self.distribute_to_vertices_and_edges()
        self.update_boundary()
        self.update_extrema()



        if self.checkpoint is True:
            self.goto_latest_checkpoint()

        if skip_initial_step is False:
            yield( self.get_time() )



        while True:
            initial_time = self.get_time()

            # Cotesting point
            if self.cotesting:
                initial_time_sc = sc.get_time()
                if not numpy.allclose(initial_time, initial_time_sc):
                    raise Exception(" Error: unequal initial_time")

                
            # Cotesting point
            if self.get_timestepping_method() == 'euler':
                self.evolve_one_euler_step(yieldstep, self.finaltime)

            elif self.get_timestepping_method() == 'rk2':
                self.evolve_one_rk2_step(yieldstep, self.finaltime)

            elif self.get_timestepping_method() == 'rk3':
                self.evolve_one_rk3_step(yieldstep, self.finaltime)



            # Cotesting point
            self.apply_fractional_steps()

            self.set_time(initial_time + self.timestep)
            # Cotesting point
            if self.cotesting:
                if not numpy.allclose(self.timestep, sc.timestep):
                    raise Exception("Error: unequal timestep %lf,  %lf" % \
                            (self.timestep, sc.timestep))
                sc.set_time(initial_time_sc + sc.timestep)



            # Cotesting point
            self.update_ghosts()

            # Cotesting point
            self.distribute_to_vertices_and_edges()

            # Cotesting point
            self.update_boundary()

            # Cotesting point
            self.update_other_quantities()

            # Cotesting point
            self.update_extrema()         



            self.number_of_steps += 1
            if self._order_ == 1:
                self.number_of_first_order_steps += 1

            if self.finaltime is not None and \
                        self.get_time() >= self.finaltime-epsilon:
                
                if self.get_time() > self.finaltime:
                    raise Exception()

                self.set_time(self.finaltime)
                self.log_operator_timestepping_statistics()
                yield(self.get_time())
                break

            if self.get_time() >= self.yieldtime:
                if self.checkpoint is True:
                    self.store_checkpoint()
                    self.delete_old_checkpoints()

                self.log_operator_timestepping_statistics()

                # Cotesting point
                if self.cotesting:
                    sc.number_of_steps += 1
                    if sc._order_ == 1:
                        sc.number_of_first_order_steps += 1

                    if sc.finaltime is not None and \
                                sc.get_time() >= sc.finaltime-epsilon:
                        if sc.get_time() > sc.finaltime:
                            raise Exception()
                        sc.set_time(sc.finaltime)
                        sc.log_operator_timestepping_statistics()
                    if sc.get_time() >= sc.yieldtime:
                        if sc.checkpoint is True:
                            sc.store_checkpoint()
                            sc.delete_old_checkpoints()

                        sc.log_operator_timestepping_statistics()


                yield(self.get_time())


                self.yieldtime += yieldstep 
                self.recorded_min_timestep = self.evolve_max_timestep
                self.recorded_max_timestep = self.evolve_min_timestep
                self.number_of_steps = 0
                self.number_of_first_order_steps = 0


                if self.using_gpu:
                    drv.memset_d32( self.max_speed_gpu, 
                                0, self.number_of_elements*2)
                else:
                    self.max_speed = numpy.zeros(N, numpy.float)




                # Cotesting point
                if self.cotesting:
                    sc.yieldtime += yieldstep 
                    sc.recorded_min_timestep = sc.evolve_max_timestep
                    sc.recorded_max_timestep = sc.evolve_min_timestep
                    sc.number_of_steps = 0
                    sc.number_of_first_order_steps = 0
                    sc.max_speed = numpy.zeros(N, numpy.float)


        ctx.synchronize()

        fin_evo = time.time()

        # FIXME: copy back
        fin_time = time.time()

        print "\nData copy in time  %lf\n" % (ini_evo - ini_time)
        print "Evolve time        %lf\n" % (fin_evo - ini_evo)
        print "Data copy out time %lf\n" % (fin_time - fin_evo)
        print "Whole time:        %lf\n" % (fin_time - ini_time)
        """ Pop up stack memory
        
        If not using PyCUDA auto_init context, we need to pop up 
        stack memory manually
        """

        global auto_init_context
        if not auto_init_context:
            ctx.pop()



    
""" Below are some utilities
"""

def get_sourceModule(k_dir, k_name, rearranged_domain=False):
    if rearranged_domain:
        defince_macro = "#define REARRANGED_DOMAIN\n"
    else:
        defince_macro = ""
    return SourceModule(
            defince_macro + open( k_dir + k_name, "r").read(),
            arch = 'compute_20',
            code = 'sm_20',
            options =['-use_fast_math', '--compiler-options', '-O3'],
            include_dirs=[ k_dir ]
            )



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



def get_device_array(a, trans=False):
    tmp_dev_p = drv.mem_alloc(a.nbytes)
    
    if trans:
        asy_cpy(a, tmp_dev_p)
    return tmp_dev_p


# FIXME:
def asy_cpy(a, a_gpu, async=False):
    global auto_init_context

    if auto_init_context and async:
        strm = drv.Stream()
        drv.memcpy_htod_async(a_gpu, a, strm)
        
        # Test correctness
        #ctx.synchronize()
        #b= numpy.zeros_like(a, a.dtype)
        #drv.memcpy_dtoh(b, a_gpu)
        #print numpy.allclose(a, b)
        return strm
    else:
        drv.memcpy_htod(a_gpu, a)



def cpy_back(a, a_gpu, async=True):
    global auto_init_context

    if auto_init_context and async:
        strm = drv.Stream()
        drv.memcpy_dtoh_async(a, a_gpu, strm)
        return strm
    else:
        drv.memcpy_dtoh(a, a_gpu)



def cpy_back_and_cmp(a, b, value_type, gpu = True, rg = False):
    if False:
        if value_type is "centroid_values":
            cpy_back(a.centroid_values, a.centroid_values_gpu)
            return numpy.allclose(a.centroid_values, b.centroid_values)
        elif value_type is "vertex_values":
            cpy_back(a.vertex_values, a.vertex_values_gpu)
            if rg:
                return check_rearranged_array(
                        b.vertex_values, a.vertex_values, 3)
            return numpy.allclose(a.vertex_values, b.vertex_values)
        elif value_type is "boundary_values":
            cpy_back(a.boundary_values, a.boundary_values_gpu)
            return numpy.allclose(a.boundary_values, b.boundary_values)
        elif value_type is "edge_values":
            cpy_back(a.edge_values, a.edge_values_gpu)
            if rg:
                return check_rearranged_array(
                        b.edge_values, a.edge_values, 3)
            return numpy.allclose(a.edge_values, b.edge_values)
        elif value_type is "x_gradient_values":
            cpy_back(a.x_gradient, a.x_gradient_gpu)
            return numpy.allclose(a.x_gradient, b.x_gradient)
        elif value_type is "y_gradient_values":
            cpy_back(a.y_gradient, a.y_gradient_gpu)
            return numpy.allclose(a.y_gradient, b.y_gradient)
        elif value_type is "explicit_update":
            cpy_back(a.explicit_update, a.explicit_update_gpu)
            return numpy.allclose(a.explicit_update, b.explicit_update)
        elif value_type is "semi_implicit_update":
            cpy_back(a.semi_implicit_update, a.semi_implicit_update_gpu)
            return numpy.allclose(a.semi_implicit_update,
                b.semi_implicit_update)
        elif value_type is "areas":
            cpy_back(a.areas, a.areas_gpu)
            return numpy.allclose(a.areas, b.areas)
        elif value_type is "surrogate_neighbours":
            cpy_back(a.surrogate_neighbours, a.surrogate_neighbours_gpu)
            if rg:
                return check_rearranged_array(
                        b.surrogate_neighbours, a.surrogate_neighbours, 3)
            return numpy.allclose(a.surrogate_neighbours, b.surrogate_neighbours)
        elif value_type is "number_of_boundaries":
            cpy_back(a.number_of_boundaries, a.number_of_boundaries_gpu)
            return numpy.allclose(a.number_of_boundaries, 
                    b.number_of_boundaries)
        elif value_type is "centroid_coordinates":
            cpy_back(a.centroid_coordinates, a.centroid_coordinates_gpu)
            if rg:
                return check_rearranged_array(
                        b.centroid_coordinates, a.centroid_coordinates, 2)
            return numpy.allclose(a.centroid_coordinates, 
                    b.centroid_coordinates)
        elif value_type is "vertex_coordinates":
            cpy_back(a.vertex_coordinates, a.vertex_coordinates_gpu)
            if rg:
                return check_rearranged_array(
                        b.vertex_coordinates, a.vertex_coordinates, 32)
            return numpy.allclose(a.vertex_coordinates, 
                    b.vertex_coordinates)
        elif value_type is "edge_coordinates":
            cpy_back(a.edge_coordinates, a.edge_coordinates_gpu)
            if rg:
                return check_rearranged_array(
                        b.edge_coordinates, a.edge_coordinates, 32)
            return numpy.allclose(a.edge_coordinates, 
                    b.edge_coordinates)
        else:
            raise Exception('Unknown value_type %s' % value_type)
    else:
        if value_type is "centroid_values":
            return numpy.allclose(a.centroid_values, b.centroid_values)
        elif value_type is "vertex_values":
            return numpy.allclose(a.vertex_values, b.vertex_values)
        elif value_type is "boundary_values":
            return numpy.allclose(a.boundary_values, b.boundary_values)
        elif value_type is "edge_values":
            return numpy.allclose(a.edge_values, b.edge_values)
        elif value_type is "x_gradient_values":
            return numpy.allclose(a.x_gradient, b.x_gradient)
        elif value_type is "y_gradient_values":
            return numpy.allclose(a.y_gradient, b.y_gradient)
        elif value_type is "explicit_update":
            return numpy.allclose(a.explicit_update, b.explicit_update)
        elif value_type is "semi_implicit_update":
            return numpy.allclose(
                a.semi_implicit_update, b.semi_implicit_update)
        elif value_type is "vertex_coordinates":
            return numpy.allclose(
                a.vertex_coordinates, b.vertex_coordinates)
        elif value_type is "areas":
            return numpy.allclose(a.areas, b.areas)
        elif value_type is "surrogate_neighbours":
            return numpy.allclose(
                a.surrogate_neighbours, b.surrogate_neighbours)
        elif value_type is "number_of_boundaries":
            return numpy.allclose(
                a.number_of_boundaries, b.number_of_boundaries)
        elif value_type is "centroid_coordinates":
            return numpy.allclose(
                a.centroid_coordinates, b.centroid_coordinates)
        elif value_type is "vertex_coordinates":
            return numpy.allclose(
                a.vertex_coordinates, b.vertex_coordinates)
        else:
            raise Exception('Unknown value_type %s' % value_type)
        




""" Below are test functions for each main cotesting level
"""

def test_distribute_to_vertexs_and_edges(domain, IO = 'Output'):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain

    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']

    s2 = domain.cotesting_domain.quantities['stage']
    xm2 = domain.cotesting_domain.quantities['xmomentum']
    ym2 = domain.cotesting_domain.quantities['ymomentum']
    e2 = domain.cotesting_domain.quantities['elevation']

    h1 = domain.quantities['height']
    xv1 = domain.quantities['xvelocity']
    yv1 = domain.quantities['yvelocity']

    h2 = sc.quantities['height']
    xv2 = sc.quantities['xvelocity']
    yv2 = sc.quantities['yvelocity']
    

    res = []


    res.append( numpy.allclose( domain.flux_timestep,
        domain.cotesting_domain.flux_timestep))
    res.append( cpy_back_and_cmp( s1, s2, 'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'x_gradient_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'y_gradient_values', gpu, rg))

    res.append( cpy_back_and_cmp( xm1, xm2,'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'x_gradient_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'y_gradient_values', gpu, rg))

    res.append( cpy_back_and_cmp( ym1, ym2,'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'x_gradient_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'y_gradient_values', gpu, rg))

    res.append( cpy_back_and_cmp( e1, e2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( e1, e2, 'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp(h1, h2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp(xv1, xv2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp(yv1, yv2,'vertex_values', gpu, rg))

    if res.count(True) != res.__len__():
        raise Exception( " --> distribute_to_vertices_and_edges ", res)



def test_evolve_one_euler_step(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain

    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']
    h1 = domain.quantities['height']
    xv1 = domain.quantities['xvelocity']
    yv1 = domain.quantities['yvelocity']
    f1 = domain.quantities['friction']

    s2 = domain.cotesting_domain.quantities['stage']
    xm2 = domain.cotesting_domain.quantities['xmomentum']
    ym2 = domain.cotesting_domain.quantities['ymomentum']
    e2 = domain.cotesting_domain.quantities['elevation']
    h2 = domain.cotesting_domain.quantities['height']
    xv2 = domain.cotesting_domain.quantities['xvelocity']
    yv2 = domain.cotesting_domain.quantities['yvelocity']
    f2 = domain.cotesting_domain.quantities['friction']


    res = []
    res.append( numpy.allclose( domain.flux_timestep,
        domain.cotesting_domain.flux_timestep))
    res.append( numpy.allclose( domain.recorded_max_timestep,
        domain.cotesting_domain.recorded_max_timestep))
    res.append( numpy.allclose( domain.recorded_min_timestep,
        domain.cotesting_domain.recorded_min_timestep))
    
    res.append( numpy.allclose( domain.smallsteps, 
        domain.cotesting_domain.smallsteps ))

    res.append( cpy_back_and_cmp( s1, s2, 'explicit_update' , gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'semi_implicit_update', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values' , gpu, rg))


    res.append( cpy_back_and_cmp( xm1, xm2,'explicit_update', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'semi_implicit_update', gpu,rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu, rg))
    

    res.append( cpy_back_and_cmp( ym1, ym2,'explicit_update', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'semi_implicit_update', gpu,rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu, rg))
    

    if res.count(True) != res.__len__():
        raise Exception( " --> evolve_one_euler_step ", res)
    
    

def test_update_ghosts(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain

    sc = domain.cotesting_domain

    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']

    s2 = domain.cotesting_domain.quantities['stage']
    xm2 = domain.cotesting_domain.quantities['xmomentum']
    ym2 = domain.cotesting_domain.quantities['ymomentum']
    e2 = domain.cotesting_domain.quantities['elevation']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values' , gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu, rg))

    # This for update_timestep check point
    res.append( cpy_back_and_cmp( e1, e2, 'edge_values' , gpu, rg))


    if res.count(True) != res.__len__():
        print " --> update_ghosts ",res



def test_update_extrema(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain

    res = []
    if res.count(True) != res.__len__():
        raise Exception( " --> update_extrema ", res)
    


def test_update_timestep(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain

    res = []
    res.append( numpy.allclose(domain._order_ , sc._order_ ))
    res.append( numpy.allclose(domain.default_order , sc.default_order ))
    res.append( numpy.allclose(domain.get_time() , sc.get_time() ))
    res.append( numpy.allclose(domain.CFL , sc.CFL ))
    res.append( numpy.allclose(domain.smallsteps , sc.smallsteps ))
    res.append( numpy.allclose(domain.max_smallsteps , sc.max_smallsteps ))
    res.append( numpy.allclose(
        domain.recorded_max_timestep , sc.recorded_max_timestep ))
    res.append( numpy.allclose(
        domain.recorded_min_timestep , sc.recorded_min_timestep ))
    res.append( numpy.allclose(
        domain.evolve_max_timestep , sc.evolve_max_timestep ))
    res.append( numpy.allclose(
        domain.evolve_min_timestep , sc.evolve_min_timestep ))
    res.append( numpy.allclose(domain.flux_timestep , sc.flux_timestep ))


    if res.count(True) != res.__len__():
        raise Exception( " --> update_timestep ",res)



def test_update_conserved_quantities(domain, output =True):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain

    for name in domain.conserved_quantities:
        Q1 = domain.quantities[name]
    
        Q2 = domain.cotesting_domain.quantities[name]
    
        res = []
        res.append( cpy_back_and_cmp( Q1, Q2, "centroid_values", gpu, rg))
        if output:
            res.append( cpy_back_and_cmp( 
                Q1, Q2, "semi_implicit_update", gpu, rg))
        res.append( cpy_back_and_cmp( Q1, Q2, "explicit_update", gpu, rg))
        res.append( numpy.allclose(
                domain.timestep, domain.cotesting_domain.timestep))
        
        if res.count(True)  != res.__len__():
            if not res[0]:
                cnt = 0
                for i in range(Q1.centroid_values.shape[0]):
                    if not numpy.allclose( Q1.centroid_values[i] , 
                                            Q2.centroid_values[i]):
                        if cnt < 5 :
                            print i, Q1.centroid_values[i], \
                                Q2.centroid_values[i]
                        cnt += 1
                print 0, cnt, Q1.centroid_values, Q2.centroid_values
            if not res[1]:
                print 1, Q1.semi_implicit_update, Q2.semi_implicit_update
            if not res[2]:
                print 2, Q1.explicit_update, Q2.explicit_update

            raise Exception("Error: update_conserved_quantities",name, res)



def test_manning_friction_implicit(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain


    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']
    f1 = domain.quantities['friction']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']
    e2 = sc.quantities['elevation']
    f2 = sc.quantities['friction']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu, rg))

    res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'semi_implicit_update', gpu,rg))

    res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'semi_implicit_update', gpu,rg))

    res.append( cpy_back_and_cmp( e1, e2, 'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp( f1, f2, 'centroid_values', gpu, rg))

    if res.count(True) != res.__len__():
        import math
        h = []
        S = []
        tmp_xmom = []
        tmp_ymom = []
        for i in range(domain.number_of_elements):
            if f2.centroid_values[i] > domain.minimum_allowed_height :
                h.append( s2.centroid_values[i] - \
                    (   e2.vertex_values[i][0] + \
                        e2.vertex_values[i][1] + \
                        e2.vertex_values[i][2])/3
                    )
                if h[i] >= domain.minimum_allowed_height: 
                    S.append( -domain.g * \
                        f1.centroid_values[i] *f1.centroid_values[i] * \
                        math.sqrt( 
                            xm2.centroid_values[i]*xm2.centroid_values[i]+\
                            ym2.centroid_values[i]*ym2.centroid_values[i]))
                    S[i] /= pow (h[i], 7.0/3)
                else:
                    S.append( 0 )
            else:
                h.append(0)
                S.append(0)

            
        if domain.use_sloped_mannings:
            raise Exception( " --> manning friction sloped ", res)
        else:
            raise Exception( " --> manning friction flat ", res)



def test_update_boundary(domain, inputOnly=False):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain

    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']
    h1 = domain.quantities['height']
    xv1 = domain.quantities['xvelocity']
    yv1 = domain.quantities['yvelocity']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']
    e2 = sc.quantities['elevation']
    h2 = sc.quantities['height']
    xv2 = sc.quantities['xvelocity']
    yv2 = sc.quantities['yvelocity']
    
    for tag in domain.tag_boundary_cells:
        B1 = domain.boundary_map[tag]
        if B1 is None :
            continue

        ids1 = domain.tag_boundary_cells[tag]
        if ids1 is None:
            continue

        B2 = sc.boundary_map[tag]

        if isinstance( B1, Reflective_boundary):
            
            ipt = []
            res = []

            ipt.append( cpy_back_and_cmp(s1, s2,'edge_values', gpu, rg))
            ipt.append( cpy_back_and_cmp(e1, e2,'edge_values', gpu, rg))
            ipt.append( cpy_back_and_cmp(h1, h2,'edge_values', gpu, rg))
            ipt.append( cpy_back_and_cmp(xm1, xm2,'edge_values', gpu, rg))
            ipt.append( cpy_back_and_cmp(ym1, ym2,'edge_values', gpu, rg))
            ipt.append( cpy_back_and_cmp(xv1, xv2,'edge_values', gpu, rg))
            ipt.append( cpy_back_and_cmp(yv1, yv2,'edge_values', gpu, rg))

            res.append( cpy_back_and_cmp(s1, s2,'boundary_values', gpu,rg))
            res.append( cpy_back_and_cmp(e1, e2,'boundary_values', gpu,rg))
            res.append( cpy_back_and_cmp(h1, h2,'boundary_values', gpu,rg))
            res.append( cpy_back_and_cmp(xm1,xm2,'boundary_values',gpu,rg))
            res.append( cpy_back_and_cmp(ym1,ym2,'boundary_values',gpu,rg))
            res.append( cpy_back_and_cmp(xv1,xv2,'boundary_values',gpu,rg))
            res.append( cpy_back_and_cmp(yv1,yv2,'boundary_values',gpu,rg))

            #print " boundary_values from refl %s" % s1.boundary_values
            #print "     -> %s " % s2.boundary_values
            #b = numpy.zeros_like(ids1)
            #cpy_back(b, domain.boundary_index[tag][0] )
            #print "   ids %s" % b, ids1
            #b = numpy.zeros_like(sc.boundary_cells[ids1])
            #cpy_back(b, domain.boundary_index[tag][1] )
            #print "   vol %s" % b, sc.boundary_cells[ids1]
            #b = numpy.zeros_like(sc.boundary_edges[ids1])
            #cpy_back(b, domain.boundary_index[tag][2] )
            #print "   edgd %s" % b, sc.boundary_edges[ids1]
            #print "   P %s" % (sc.boundary_cells[ids1] *3 + sc.boundary_edges[ids1])
            if ipt.count(True) != ipt.__len__():
               print "\n  Input values check ", ipt
               if not res[0]:
                   print "se", s1.edge_values, s2.edge_values
               if not res[1]:
                   print "ee", e1.edge_values, e2.edge_values
               if not res[2]:
                   print "he", h1.edge_values, h2.edge_values
               if not res[3]:
                   print "xme", xm1.edge_values, xm2.edge_values
               if not res[4]:
                   print "yme", ym1.edge_values, ym2.edge_values
               if not res[5]:
                   print "xve", xv1.edge_values, xv2.edge_values
               if not res[6]:
                   print "yve", yv1.edge_values, yv2.edge_values

            if not inputOnly and res.count(True) != res.__len__():
                print "\n Error in update_boundary", tag, res
                if not res[0]:
                    print "sb", s1.boundary_values, s2.boundary_values
                if not res[1]:
                    print "eb", e1.boundary_values, e2.boundary_values
                if not res[2]:
                    print "hb", h1.boundary_values, h2.boundary_values
                if not res[3]:
                    print "xmb", xm1.boundary_values, xm2.boundary_values
                if not res[4]:
                    print "ymb", ym1.boundary_values, ym2.boundary_values
                if not res[5]:
                    print "xvb", xv1.boundary_values, xv2.boundary_values
                if not res[6]:
                    print "yvb", yv1.boundary_values, yv2.boundary_values



                raise Exception("Error in update_boundary reflective "+ tag)

        elif isinstance( B1, Dirichlet_boundary ):
            for j ,name in enumerate(domain.evolved_quantities):
                Q1 = domain.quantities[name]
                Q2 = sc.quantities[name]
                if not cpy_back_and_cmp(Q1, Q2, "boundary_values"):
                    print "Error in update_boundary", tag, name, Q1.boundary_values, Q2.boundary_values
                    #cpy_back(Q1.edge_values, Q1.edge_values_gpu)
                    #print domain.quantities[name].edge_values
                    #print sc.quantities[name].edge_values 
                    raise Exception()
            if len( B1.dirichlet_values ) != \
                        len(domain.evolved_quantities):
                    
                for j, name in enumerate(domain.conserved_quantities):
                    Q1 = domain.quantities[name].boundary_values
                    Q2 = sc.quantities[name].boundary_values
                    if not numpy.allclose(Q1, Q2):
                        print "Error in update_boundary", tag, name
                        raise Exception()


def test_update_other_quantities(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain

    h1 = domain.quantities['height']
    xv1 = domain.quantities['xvelocity']
    yv1 = domain.quantities['yvelocity']

    h2 = sc.quantities['height']
    xv2 = sc.quantities['xvelocity']
    yv2 = sc.quantities['yvelocity']
    
    res = []

    res.append( cpy_back_and_cmp(h1, h2,'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp(xv1, xv2,'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp(yv1, yv2,'edge_values', gpu, rg))

    res.append( cpy_back_and_cmp(h1, h2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp(xv1, xv2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp(yv1, yv2,'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp(h1, h2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp(xv1, xv2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp(yv1, yv2,'centroid_values', gpu, rg))

    if res.count(True) != res.__len__():
       raise Exception("Error in update_other_quantities", res)


def test_compute_fluxes(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain

    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']
    h1 = domain.quantities['height']
    xv1 = domain.quantities['xvelocity']
    yv1 = domain.quantities['yvelocity']
    f1 = domain.quantities['friction']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']
    e2 = sc.quantities['elevation']
    h2 = sc.quantities['height']
    xv2 = sc.quantities['xvelocity']
    yv2 = sc.quantities['yvelocity']
    f2 = sc.quantities['friction']


    res = []
    res.append( numpy.allclose(domain.flux_timestep, sc.flux_timestep))


    res.append( cpy_back_and_cmp( s1, s2, 'explicit_update' , gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'edge_values' , gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'boundary_values' , gpu, rg))


    res.append( cpy_back_and_cmp( xm1, xm2,'explicit_update', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'boundary_values', gpu, rg))
    

    res.append( cpy_back_and_cmp( ym1, ym2,'explicit_update', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'boundary_values', gpu, rg))
    

    if res.count(True) != res.__len__():
        if not res[0]:
            print "   flux_timestep %.9lf %.9lf" % (
                    domain.flux_timestep, sc.flux_timestep)
        if not res[1]:
            for i in range(domain.number_of_elements):
                if not numpy.allclose( 
                        s1.explicit_update[i], s2.explicit_update[i]):
                    print 0,i, s1.explicit_update[i], s2.explicit_update[i]
        if not res[4]:
            for i in range(domain.number_of_elements):
                if not numpy.allclose( 
                        xm1.explicit_update[i], xm2.explicit_update[i]):
                    print 1,i,xm1.explicit_update[i],xm2.explicit_update[i]
        if not res[7]:
            for i in range(domain.number_of_elements):
                if not numpy.allclose( 
                        ym1.explicit_update[i], ym2.explicit_update[i]):
                    print 2,i,ym1.explicit_update[i],ym2.explicit_update[i]

        res_name = ['flux_timestep',
                    'stage_explicit', 'stage_edge', 'stage_boundary', 
                    'xmom_explicit', 'xmom_edge', 'xmom_boundary',
                    'ymom_explicit', 'ymom_edge', 'ymom_boundary']
        raise Exception( " --> compute_fluxes", [ b for a,b in zip(res, res_name) if not a ])


def test_compute_forcing_terms(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain


    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']
    f1 = domain.quantities['friction']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']
    e2 = sc.quantities['elevation']
    f2 = sc.quantities['friction']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu, rg))

    res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'semi_implicit_update', gpu,rg))

    res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'semi_implicit_update', gpu,rg))

    res.append( cpy_back_and_cmp( e1, e2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( e1, e2, 'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp( f1, f2, 'centroid_values', gpu, rg))

    if res.count(True) != res.__len__():
        raise Exception( " --> compute_forcing_terms ", res)


def test_protect_against_infinitesimal_and_negative_heights(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain


    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']
    e2 = sc.quantities['elevation']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu, rg))

    res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu, rg))

    res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu, rg))

    res.append( cpy_back_and_cmp( e1, e2, 'centroid_values', gpu, rg))


    if res.count(True) != res.__len__():
        raise Exception( " --> protect_against_infinitesimal_and_negative_heights ",res)


def test_extrapolate_second_order_sw(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain


    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']
    e2 = sc.quantities['elevation']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp( e1, e2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( e1, e2, 'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'vertex_values', gpu, rg))

    res.append( domain.epsilon == sc.epsilon)
    res.append( domain.minimum_allowed_height == sc.minimum_allowed_height)
    res.append( domain.beta_w == sc.beta_w)
    res.append( domain.beta_w_dry == sc.beta_w_dry)
    res.append( domain.beta_uh == sc.beta_uh)
    res.append( domain.beta_uh_dry == sc.beta_uh_dry)
    res.append( domain.beta_vh == sc.beta_vh)
    res.append( domain.beta_vh_dry == sc.beta_vh_dry)
    res.append( domain.optimise_dry_cells == sc.optimise_dry_cells)

    res.append( cpy_back_and_cmp( domain,sc,"surrogate_neighbours",gpu,rg))
    res.append( cpy_back_and_cmp( domain,sc,"number_of_boundaries",gpu,rg))
    res.append( cpy_back_and_cmp( domain,sc,"centroid_coordinates",gpu,rg))
    res.append( cpy_back_and_cmp( domain,sc,"vertex_coordinates",gpu,rg))

    if res.count(True) != res.__len__():
        print res

        #from anuga_cuda.extrapolate.extrapolate_second_order_sw import \
        #    extrapolate_second_order_sw_python as extra 
        #

        #ss = numpy.zeros_like(s1.centroid_values, dtype=numpy.float64)
        #xs = numpy.zeros_like(s1.centroid_values, dtype=numpy.float64)
        #ys = numpy.zeros_like(s1.centroid_values, dtype=numpy.float64)

        #drv.memcpy_dtoh(ss, domain.stage_centroid_store_gpu)
        #drv.memcpy_dtoh(xs, domain.xmomentum_centroid_store_gpu)
        #drv.memcpy_dtoh(ys, domain.ymomentum_centroid_store_gpu)

        #extra(sc, ss, xs, ys)


        #gpu = False
        #res = []
        #res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu))
        #res.append( cpy_back_and_cmp( s1, s2, 'vertex_values', gpu))

        #res.append( cpy_back_and_cmp( e1, e2, 'centroid_values', gpu))
        #res.append( cpy_back_and_cmp( e1, e2, 'vertex_values', gpu))

        #res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu))
        #res.append( cpy_back_and_cmp( xm1, xm2,'vertex_values', gpu))

        #res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu))
        #res.append( cpy_back_and_cmp( ym1, ym2,'vertex_values', gpu))

        #res.append( cpy_back_and_cmp( domain, sc, 
        #        "surrogate_neighbours", gpu))
        #res.append( cpy_back_and_cmp( domain, sc, 
        #        "number_of_boundaries", gpu))
        #res.append( cpy_back_and_cmp( domain, sc, 
        #        "centroid_coordinates", gpu))
        #res.append( cpy_back_and_cmp( domain, sc, 
        #        "vertex_coordinates", gpu))

        #print res
        #cnt = 0
        #if not res[5]:
        if False:
            for i in range(domain.number_of_elements):
                if (xm1.vertex_values[i] != xm2.vertex_values[i]).all():
                    if domain.number_of_boundaries[i] == 1:
                        print i, xm1.vertex_values[i], xm2.vertex_values[i]
                    cnt += 1
            print cnt

        raise Exception( " --> extrapolate_second_order_sw ", res)


def test_balance_deep_and_shallow(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain


    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']
    e2 = sc.quantities['elevation']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp( xm1, xm2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp( ym1, ym2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'vertex_values', gpu, rg))

    res.append( cpy_back_and_cmp( e1, e2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( e1, e2, 'vertex_values', gpu, rg))


    if res.count(True) != res.__len__():
        raise Exception( " --> balance_deep_and_shallow ", res)
        

def test_interpolate_from_vertices_to_edges(domain):
    # For time-dependence issues
    ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain


    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'edge_values', gpu, rg))

    res.append( cpy_back_and_cmp( xm1, xm2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'edge_values', gpu, rg))

    res.append( cpy_back_and_cmp( ym1, ym2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'edge_values', gpu, rg))


    if res.count(True) != res.__len__():
        for i in range(domain.number_of_elements):
            if not numpy.allclose(s1.edge_values[i], s2.edge_values[i]):
                print i, s1.edge_values[i], s2.edge_values[i]
        raise Exception( " --> interpolate_from_vertices_to_edges ", res)


def test_extrapolate_second_order_and_limit_by_vertex(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain


    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'x_gradient', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'y_gradient', gpu, rg))


    res.append( cpy_back_and_cmp( x1, x2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( x1, x2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( x1, x2, 'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp( x1, x2, 'x_gradient', gpu, rg))
    res.append( cpy_back_and_cmp( x1, x2, 'y_gradient', gpu, rg))


    res.append( cpy_back_and_cmp( y1, y2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( y1, y2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( y1, y2, 'edge_values', gpu, rg))
    res.append( cpy_back_and_cmp( y1, y2, 'x_gradient', gpu, rg))
    res.append( cpy_back_and_cmp( y1, y2, 'y_gradient', gpu, rg))


    if res.count(True) != res.__len__():
        raise Exception( " --> extrapolate_second_order_and_limit_by_vertex ", res)

    
def test_extrapolate_first_order(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain


    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']

    res = []
    res.append( cpy_back_and_cmp( s1, s2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'edge_values', gpu, rg))


    res.append( cpy_back_and_cmp( x1, x2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( x1, x2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( x1, x2, 'edge_values', gpu, rg))

    res.append( cpy_back_and_cmp( y1, y2, 'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp( y1, y2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( y1, y2, 'edge_values', gpu, rg))


    if res.count(True) != res.__len__():
        raise Exception( " --> extrapolate_first_order ", res)


def test_update_centroids_of_velocities_and_height(domain):
    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain

    s1 = domain.quantities['stage']
    xm1 = domain.quantities['xmomentum']
    ym1 = domain.quantities['ymomentum']
    e1 = domain.quantities['elevation']
    h1 = domain.quantities['height']
    xv1 = domain.quantities['xvelocity']
    yv1 = domain.quantities['yvelocity']

    s2 = sc.quantities['stage']
    xm2 = sc.quantities['xmomentum']
    ym2 = sc.quantities['ymomentum']
    e2 = sc.quantities['elevation']
    h2 = sc.quantities['height']
    xv2 = sc.quantities['xvelocity']
    yv2 = sc.quantities['yvelocity']
    
            
    res = []

    res.append( cpy_back_and_cmp(s1, s2,'boundary_values', gpu, rg))
    res.append( cpy_back_and_cmp(e1, e2,'boundary_values', gpu, rg))
    res.append( cpy_back_and_cmp(h1, h2,'boundary_values', gpu, rg))
    res.append( cpy_back_and_cmp(xm1, xm2,'boundary_values', gpu, rg))
    res.append( cpy_back_and_cmp(ym1, ym2,'boundary_values', gpu, rg))
    res.append( cpy_back_and_cmp(xv1, xv2,'boundary_values', gpu, rg))
    res.append( cpy_back_and_cmp(yv1, yv2,'boundary_values', gpu, rg))

    res.append( cpy_back_and_cmp(s1, s2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp(e1, e2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp(h1, h2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp(xm1, xm2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp(ym1, ym2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp(xv1, xv2,'centroid_values', gpu, rg))
    res.append( cpy_back_and_cmp(yv1, yv2,'centroid_values', gpu, rg))

    if res.count(True) != res.__len__():
        if not res[12]:
            for i in range(domain.number_of_elements):
                if not numpy.allclose(
                        xv1.centroid_values[i],
                        xv2.centroid_values[i]) :
                    print "  xvelocity centroid %d %lf %lf" % \
                        (i, xv1.centroid_values[i], xv2.centroid_values[i])
        raise Exception("Error in update_centroids_of_velocities_and_height", res)


