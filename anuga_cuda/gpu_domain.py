#!/usr/bin/env python

# System module
import numpy 

# ANUGA module
from anuga import Domain
from anuga.abstract_2d_finite_volumes.generic_domain \
                    import Generic_Domain

from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import \
        Transmissive_boundary, Dirichlet_boundary, Compute_fluxes_boundary, \
        Time_boundary, Time_boundary, Time_boundary
from anuga.shallow_water.boundaries import Reflective_boundary

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


# Config info
#from anuga_cuda.config import compute_fluxes_dir
#from anuga_cuda.config import gravity_dir
#from anuga_cuda.config import extrapolate_dir
#from anuga_cuda.config import protect_dir
#from anuga_cuda.config import balance_dir
#from anuga_cuda.config import interpolate_dir
#from anuga_cuda.config import evaluate_dir
#from anuga_cuda.config import get_absolute_dir
#from anuga_cuda.config import manning_friction_dir
#from anuga_cuda.config import saxpy_dir
#from anuga_cuda.config import set_boundary_dir
#from anuga_cuda.config import update_centroids_dir
from anuga_cuda import kernel_path as kp



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
            using_gpu=False): 

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

        self.boundary_index = {}


        self.using_gpu = using_gpu

        #self.end_event = drv.Event()

        if using_page_locked:
            self.timestep_array = drv.pagelocked_zeros(
                    self.number_of_elements,
                    dtype = numpy.float64,
                    mem_flags=drv.host_alloc_flags.DEVICEMAP)
        else:
            self.timestep_array = numpy.zeros(
                    self.number_of_elements,
                    dtype = numpy.float64)
       
        # compute_fluxes function
        self.compute_fluxes_mod = SourceModule(
                open( kp["compute_fluxes_dir"] + \
                        "compute_fluxes.cu","r" ).read(),
                include_dirs=[kp["compute_fluxes_dir"]]
                )

        self.compute_fluxes_func = self.compute_fluxes_mod.get_function(
                #"compute_fluxes_central_structure_cuda_single")
                "compute_fluxes_central_structure_CUDA")

        # gravity_wb function
        self.gravity_wb_mod = SourceModule(
                open( kp["gravity_dir"]+"gravity.cu","r ").read(),
                include_dirs=[ kp["gravity_dir"] ]
                )

        self.gravity_wb_func = self.gravity_wb_mod.get_function(
                "gravity_wb")

        # extrapolate_second_order_sw function
        self.extrapolate_second_order_sw_mod = SourceModule(
                open( kp["extrapolate_dir"] +\
                        "extrapolate_second_order_sw.cu").read(),
                options=["--ptxas-options=-v"],
                include_dirs=[ kp["extrapolate_dir"]]
            )

        self.extrapolate_second_order_sw_true_func = \
            self.extrapolate_second_order_sw_mod.get_function(
                    "extrapolate_second_order_sw_true")

        self.extrapolate_second_order_sw_false_func = \
            self.extrapolate_second_order_sw_mod.get_function(
                "extrapolate_second_order_sw_false")

        # extrapolate_second_order_and_limit_by_vertex_or_edge function
        self.extrapolate_second_order_and_limit_by_vertex_or_edge_mod = \
            SourceModule(
                open( kp["extrapolate_dir"] + \
                "extrapolate_second_order_and_limit_by_vertex_or_edge.cu"
                    ).read(),
                include_dirs=[ kp["extrapolate_dir"]]
                )

        self.extrapolate_second_order_and_limit_by_vertex_func = \
            self.extrapolate_second_order_and_limit_by_vertex_or_edge_mod.\
                get_function(
                    "extrapolate_second_order_and_limit_by_vertex")

        self.extrapolate_second_order_and_limit_by_edge_func = \
            self.extrapolate_second_order_and_limit_by_vertex_or_edge_mod.\
                get_function(
                    "extrapolate_second_order_and_limit_by_edge")

        # extrapolate_first_order function
        self.extrapolate_first_order_mod = \
            SourceModule(
                open( kp["extrapolate_dir"] +\
                    "extrapolate_first_order.cu").read(),
                include_dirs=[ kp["extrapolate_dir"] ]
                )

        self.extrapolate_first_order_func = \
            self.extrapolate_first_order_mod.get_function(
                    "extrapolate_first_order")


        # protect function
        self.protect_mod = SourceModule(
                open( kp["protect_dir"] + "protect.cu").read(),
                include_dirs=[ kp["protect_dir"] ]
                )

        self.protect_sw_func = self.protect_mod.get_function(
                "_protect_sw")
        self.protect_swb2_func = self.protect_mod.get_function(
                "_protect_swb2")
        
        # balance function
        self.balance_deep_and_shallow_mod = SourceModule(
                open( kp["balance_dir"] + "balance_deep_and_shallow.cu").read(),
                include_dirs=[ kp["balance_dir"]]
                )

        self.balance_deep_and_shallow_func = \
            self.balance_deep_and_shallow_mod.get_function(
                "_balance_deep_and_shallow")

        # interpolate_from_vertices_to_edges function
        self.interpolate_from_vertices_to_edges_mod = SourceModule(
                open( kp["interpolate_dir"] + \
                    "interpolate_from_vertices_to_edges.cu").read(),
                include_dirs=[ kp["interpolate_dir"]]
                )

        self.interpolate_from_vertices_to_edges_func = \
            self.interpolate_from_vertices_to_edges_mod.get_function(
                "_interpolate_from_vertices_to_edges")


        # evaluate_segment function
        self.evaluate_segment_mod = SourceModule(
                open( kp["evaluate_dir"] +"evaluate_segment.cu").read(),
                include_dirs=[ kp["evaluate_dir"] ]
                )

        self.evaluate_segment_reflective_func = \
            self.evaluate_segment_mod.get_function(
                "evaluate_segment_reflective")
        self.evaluate_segment_dirichlet_1_func = \
            self.evaluate_segment_mod.get_function(
                "evaluate_segment_dirichlet_1")
        self.evaluate_segment_dirichlet_2_func = \
            self.evaluate_segment_mod.get_function(
                "evaluate_segment_dirichlet_2")


        # get_absolute function
        self.get_absolute_mod = SourceModule(
                open( kp["get_absolute_dir"] +"get_absolute.cu").read(),
                include_dirs=[ kp["get_absolute_dir"] ]
                )
            
        self.get_absolute_func = \
            self.get_absolute_mod.get_function("get_absolute")
        

        # manning_friction function
        self.manning_friction_mod = SourceModule(
                open( kp["manning_friction_dir"] + \
                    "manning_friction.cu").read(),
                include_dirs=[ kp["manning_friction_dir"] ]
                )

        self.manning_friction_sloped_func = \
            self.manning_friction_mod.get_function(
                "_manning_friction_sloped")

        self.manning_friction_flat_func = \
            self.manning_friction_mod.get_function(
                "_manning_friction_flat")

        
        # saxpy_centroid_values function
        self.saxpy_centroid_values_mod = SourceModule(
                open( kp["saxpy_dir"] + "saxpy_centroid_values.cu").read(),
                include_dirs=[ kp["saxpy_dir"] ]
                )

        self.saxpy_centroid_values_func = \
            self.saxpy_centroid_values_mod.get_function(
                "_saxpy_centroid_values")

        
        # set_boundry function
        self.set_boundary_mod = SourceModule(
                open( kp["set_boundary_dir"] + "set_boundary.cu").read(),
                include_dirs=[ kp["set_boundary_dir"] ]
                )
                
        self.set_boundary_values_from_edges_func = \
            self.set_boundary_mod.get_function(
                "set_boundary_values_from_edges")

        # update_centroids_of_velocities_and_height function
        self.update_centroids_of_velocities_and_height_mod = SourceModule(
                open( kp["update_centroids_dir"] +\
                    "update_centroids_of_velocities_and_height.cu").read(),
                include_dirs=[ kp["update_centroids_dir"] ]
                )
        self.update_centroids_of_velocities_and_height_func = \
            self.update_centroids_of_velocities_and_height_mod.\
            get_function("update_centroids_of_velocities_and_height")

        
        # update function
        self.update_mod = SourceModule(
                open(kp["update_dir"]+"update.cu").read(),
                include_dirs=[ kp["update_dir"] ])
        self.update_func = self.update_mod.get_function("update")

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



    def asynchronous_transfer(self):
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


    def compute_fluxes(self):
        if self.using_gpu :
            W1 = 32
            W2 = 1
            W3 = 1
            

            #strm_cf = drv.Stream()

            self.compute_fluxes_func(
                numpy.int32(self.number_of_elements),
                numpy.float64(self.evolve_max_timestep),
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
                grid = ((self.number_of_elements+W1*W2*W3-1)/(W1*W2*W3), 1)
                #stream = strm_cf
                )
                
            #strm_g = drv.Stream()
            drv.memcpy_dtoh(self.timestep_array, self.timestep_array_gpu)

            self.gravity_wb_func(
                numpy.uint64(self.number_of_elements),
                numpy.float64(self.g),
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
                block = (W1, W2, W3),
                grid =((self.number_of_elements + W1*W2*W3-1)/(W1*W2*W3),1)
                #stream =strm_g
                )

            #drv.memcpy_dtoh_async(
            #        self.timestep_array, 
            #        self.timestep_array_gpu, 
            #        strm_cf
            #        )

            
            
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
                numpy.int32(self.number_of_elements),
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
                        self.mim_elevation_edgevalue_gpu,
                        self.max_elevation_edgevalue_gpu,
                        self.count_wet_neighbours_gpu,
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
                            numpy.int32(self.number_of_elements),
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
                        self.extrapolate_second_order_and_limit_by_edge_func(
                                numpy.float64( Q.beta),
                                self.centroid_coordinates_gpu,
                                self.vertex_coordinates_gpu,
                                self.number_of_boundaries_gpu,
                                self.surrogate_neighbours_gpu,
                                self.neighbours_gpu,

                                Q.centroid_values_gpu,
                                Q.vertex_values_gpu,
                                Q.edge_values_gpu,
                                Q.x_gradient_gpu,
                                Q.y_gradient_gpu,
                                block = (W1, W2, W3),
                                grid=((self.number_of_elements+W1*W2*W3-1)/(W1*W2),1)
                                )
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
                            #Q.extrapolate_first_order()
                            self.extrapolate_first_order_func(
                                numpy.int32(self.number_of_elements),
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
                        self.extrapolate_second_order_sw()
                    else:
                        raise Exception('Unknown order')

                else:
                    for name in self.conserved_quantities:
                        Q = self.quantities[name]

                        if self._order_ == 1:
                            #Q.extrapolate_first_order()
                            self.extrapolate_first_order_func(
                                numpy.int32(self.number_of_elements),
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
                            #Q.extrapolate_second_order_and_limit_by_vertex()
                            self.extrapolate_second_order_and_limit_by_vertex_func(
                                numpy.float64( Q.beta),
                                self.centroid_coordinates_gpu,
                                self.vertex_coordinates_gpu,
                                self.number_of_boundaries_gpu,
                                self.surrogate_neighbours_gpu,
                                self.neighbours_gpu,

                                Q.centroid_values_gpu,
                                Q.vertex_values_gpu,
                                Q.edge_values_gpu,
                                Q.x_gradient_gpu,
                                Q.y_gradient_gpu,
                                block = (W1, W2, W3),
                                grid=((self.number_of_elements+W1*W2*W3-1)/(W1*W2),1)
                                    )


                        else:
                            raise Exception('Unknown order')

                self.balance_deep_and_shallow()

                for name in self.conserved_quantities:
                    Q = self.quantities[name]
                    #Q.interpolate_from_vertices_to_edges()
                    self.interpolate_from_vertices_to_edges_func(
                            numpy.int32(Q.vertex_values.shape[0]),
                            #numpy.int32(self.number_of_elements),
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
                        numpy.int32(self.number_of_elements),
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


    def update_boundary(self):
        if self.using_gpu:
            #FIXME:result not correct
            W1 = 32
            W2 = 1
            W3 = 1

            for tag in self.tag_boundary_cells:
                B = self.boundary_map[tag]
                if B is None:
                    continue
                
                #segment_edges = self.tag_boundary_cells[tag]
                ids = self.tag_boundary_cells[tag]

                # def evaluate_segment(self, domain, segment_edges):
                if ids is None:
                    continue

                if isinstance(B, Reflective_boundary):
                    self.evaluate_segment_reflective_func(
                        numpy.int32(len(ids)),
                        self.boundary_index[tag][0],
                        self.boundary_index[tag][1],
                        self.boundary_index[tag][2],

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
                        self.quantities['yvelocity'].boundary_values_gpu,
                        
                        block = (W1, W2, W3),
                        grid=((len(ids)+W1*W2*W3-1)/(W1*W2*W3),1)
                        )

                elif isinstance(B, Dirichlet_boundary):
                    q_bdry = B.dirichlet_values
                    conserved_quantities = True
                    if len(q_bdry) == len(self.evolved_quantities):
                        conserved_quantities = False
                    if  conserved_quantities:
                        for j, name in enumerate(self.evolved_quantities):
                            Q = self.quantities[name]
                            #Q.boundary_values[ids] = Q.edge_values[vol_ids,edge_ids]
                            self.evaluate_segment_dirichlet_1_func(
                                    numpy.int32(len(ids)),
                                    self.boundary_index[tag][0],
                                    self.boundary_index[tag][1],
                                    self.boundary_index[tag][2],

                                    Q.boundary_values_gpu,
                                    Q.edge_values_gpu,
                                    block = (W1, W2, W3),
                                    grid=((len(ids)+W1*W2*W3-1)/(W1*W2*W3),1)
                                    )


                    if conserved_quantities:
                        quantities = self.conserved_quantities
                    else:
                        quantities = self.evolved_quantities

                    for j, name in enumerate(quantities):
                        Q = self.quantities[name]
                        #Q.boundary_values[ids] = q_bdry[j]
                        self.evaluate_segment_dirichlet_2_func(
                                numpy.int32(len(ids)),
                                numpy.float64(q_bdry[j]),
                                self.boundary_index[tag][0],
                                Q.boundary_values_gpu,
                                block = (W1, W2, W3),
                                grid=((len(ids)+W1*W2*W3-1)/(W1*W2*W3),1)
                                )

                else:
                    raise Exception("Can not find right type")
                
        else:
            Generic_Domain.update_boundary(self)


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
                    
                    offset=num.array([self.geo_reference.get_xllcorner(),
                        self.geo_reference.get_yllcorner()], num.float)

                    return V[i3:i3+3,:] + offset   
                else:
                    return V[i3:i3+3,:]
        else:
            return Domain.get_vertex_coordinates(self)
                    
                
    def manning_friction_implicit(self):
        """From shallow_water_domain"""  
        if self.using_gpu:
            N = self.number_of_elements
            W1 = 32
            W2 = 1
            W3 = 1
            #FIXME
            x = self.get_vertex_coordinates()
            if self.use_sloped_mannings:
                self.manning_friction_sloped_func(
                   numpy.int32(N),
                   numpy.float64(self.g),
                   numpy.float64(self.minimum_allowed_height),
   
                   self.vertex_coordinates_gpu,
                   self.quantities['stage'].centroid_values_gpu,
                   self.quantities['xmomentum'].centroid_values_gpu,
                   self.quantities['ymomentum'].centroid_values_gpu,
                   self.quantities['elevation'].vertex_values_gpu,
   
                   self.quantities['friction'].centroid_values_gpu,
                   self.quantities['xmomentum'].semi_implicit_update_gpu,
                   self.quantities['ymomentum'].semi_implicit_update_gpu,
                   block = (W1, W2, W3),
                   grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
                   )
            else:
                self.manning_friction_flat_func(
                   numpy.int32(N),
                   numpy.float64(self.g),
                   numpy.float64(self.minimum_allowed_height),
   
                   self.quantities['stage'].centroid_values_gpu,
                   self.quantities['xmomentum'].centroid_values_gpu,
                   self.quantities['ymomentum'].centroid_values_gpu,
                   self.quantities['elevation'].centroid_values_gpu,
   
                   self.quantities['friction'].centroid_values_gpu,
                   self.quantities['xmomentum'].semi_implicit_update_gpu,
                   self.quantities['ymomentum'].semi_implicit_update_gpu,
                   
                   block = (W1, W2, W3),
                   grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
                   )
   
        else:
            from anuga.shallow_water.shallow_water_domain import \
                    manning_friction_implicit
            manning_friction_implicit(self)


    def manning_friction_explicit(self):
        """From shallow_water_domain"""  
        if self.using_gpu:
           #FIXME
           x = self.get_vertex_coordinates()
   
           if self.use_sloped_mannings:
               self.manning_friction_sloped_func(
                   numpy.int32(self.number_of_elements),
                   numpy.float64(self.g),
                   numpy.float64(self.minimum_allowed_height),
   
                   self.vertex_coordinates_gpu,
                   self.quantities['stage'].centroid_values_gpu,
                   self.quantities['xmomentum'].centroid_values_gpu,
                   self.quantities['ymomentum'].centroid_values_gpu,
                   self.quantities['elevation'].vertex_values_gpu,
   
                   self.quantities['friction'].centroid_values_gpu,
                   self.quantities['xmomentum'].explicit_update_gpu,
                   self.quantities['ymomentum'].explicit_update_gpu,
                   block = (W1, W2, W3),
                   grid=((len(points.shape[0])+W1*W2*W3-1)/(W1*W2*W3),1)
                   )
           else:
               self.manning_friction_flat_func(
                   numpy.int32(self.number_of_elements),
                   numpy.float64(self.g),
                   numpy.float64(self.minimum_allowed_height),
   
                   self.quantities['stage'].centroid_values_gpu,
                   self.quantities['xmomentum'].centroid_values_gpu,
                   self.quantities['ymomentum'].centroid_values_gpu,
                   self.quantities['elevation'].vertex_values_gpu,
   
                   self.quantities['friction'].centroid_values_gpu,
                   self.quantities['xmomentum'].explicit_update_gpu,
                   self.quantities['ymomentum'].explicit_update_gpu,
                   
                   block = (W1, W2, W3),
                   grid=((len(points.shape[0])+W1*W2*W3-1)/(W1*W2*W3),1)
                   )
   
        else:
            from anuga.shallow_water.shallow_water_domain import manning_friction_implicit
            manning_friction_implicit(self)


    def compute_forcing_terms(self):
        if self.using_gpu:
            for f in self.forcing_terms:
                f()
        else:
            for f in self.forcing_terms:
                f(self)


    def update_conserved_quantities(self):
        if self.using_gpu :
            N = self.number_of_elements
            W1 = 32 
            W2 = 1
            W3 = 1
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                self.update_func(
                    numpy.float64(N),
                    numpy.float64(self.timestep),
                    Q.centroid_values_gpu,
                    Q.explicit_update_gpu,
                    Q.semi_implicit_update_gpu,
                    block = (W1, W2, W3),
                    grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
                    )

                drv.memset_d32(Q.semi_implicit_update_gpu, 0, N*2)
        else:
            Domain.update_conserved_quantities(self)


    def backup_conserved_quantities(self):
        if self.using_gpu:
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                drv.memcpy_dtod(
                        Q.centroid_backup_values_gpu, 
                        Q.centroid_values_gpu)
        else:
            Domain.backup_conserved_quantities(self)


    def saxpy_conserved_quantities(self, a, b):
        if self.using_gpu:
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                self.saxpy_centroid_values(
                    numpy.int32(self.number_of_elements),
                    numpy.float64(a),
                    numpy.float64(b),
                    Q.centroid_values_gpu,
                    Q.centroid_backup_values_gpu,
                    block = (W1, W2, W3),
                    grid=((len(points.shape[0])+W1*W2*W3-1)/(W1*W2*W3),1)
                    )
        else:
            Domain.saxpy_conserved_quantities(self, a, b)


    def update_centroids_of_velocities_and_height(self):
        if self.using_gpu:
            N = self.quantities['elevation'].boundary_values.shape[0]
            W1 = 32
            W2 = 1
            W3 = 1

            self.set_boundary_values_from_edges_func(
                numpy.int32(N),
                self.boundary_cells_gpu,
                self.boundary_edges_gpu,
                self.quantities['elevation'].boundary_values_gpu,
                self.quantities['elevation'].edge_values_gpu,
                block = (W1, W2, W3),
                grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
                )

            self.update_centroids_of_velocities_and_height_func(
                numpy.int32(N),
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
                self.quantities['yvelocity'].boundary_values_gpu,
                block = (W1, W2, W3),
                grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
                )

        else:
            Domain.update_centroids_of_velocities_and_height(self)


    def copy_back_necessary_data(self):
        pass       


    def store_timestep(self):
        if self.using_gpu:
            self.copy_back_necessary_data()
        self.writer.store_timestep()


    def evolve(self, 
                yieldstep=None,
                finaltime=None,
                duration=None,
                skip_initial_step=False):

        print " # of elements: %d" % self.number_of_elements
        if self.using_gpu:
            from anuga_cuda import sort_domain
            sort_domain(self)

            from anuga.shallow_water.shallow_water_domain import \
                    manning_friction_implicit, manning_friction_explicit

            f = self.forcing_terms
            for i in range(len(f)):
                if f[i] == manning_friction_implicit:
                    f[i] = self.manning_friction_implicit
                elif f[i] == manning_friction_explicit:
                    f[i] == self.manning_friction_explicit
                #FIXME
                # in shallow_water_domain set_gravity_methon function
                # gravity functions can be the elements of forcing_terms
                
                    
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
    global auto_init_context

    if auto_init_context:
        strm = drv.Stream()
        drv.memcpy_htod_async(a_gpu, a, strm)
        return strm
    else:
        drv.memcpy_htod(a_gpu, a)

def cpy_back(a, a_gpu):
    drv.memcpy_dtoh(a, a_gpu)
