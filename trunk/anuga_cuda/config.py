WorkingStation = "GTX480"
if WorkingStation == "Xe":
    work_dir="/home/659/zxw659/anuga-cuda/anuga_cuda/"
else:
    work_dir="/home/john/anuga_cuda/"


testing_domain_path = work_dir+"testing_domains/"

merimbula_dir = testing_domain_path+"merimbula/"

cairns_dir= testing_domain_path+"cairns/"


kernel_path = {
    "compute_fluxes_dir" : work_dir+"compute_fluxes/",
    "gravity_dir" : work_dir+"gravity/",
    "extrapolate_dir" : work_dir+"extrapolate/",
    "protect_dir" : work_dir+"protect/",
    "balance_dir" : work_dir+"balance/",
    "interpolate_dir" : work_dir+"interpolate_from_vertices_to_edges/",
    "evaluate_dir" : work_dir+"evaluate_segment/",
    "get_absolute_dir" : work_dir+"get_absolute/",
    "manning_friction_dir" : work_dir+"manning_friction/",
    "saxpy_dir" : work_dir+"saxpy_centroid_values/",
    "set_boundary_dir" : work_dir+"set_boundary/",
    "update_centroids_dir":work_dir+"update_centroids_of_velocities_and_height/",
    "update_dir" : work_dir+"update/"
    }



if WorkingStation == "Xe":
    kernel_block_configuration = {
        "balance_fun" : 128,

        "compute_fluxes_fun" : 128,

        "extrapolate_first_order_fun" : 128,

        "extrapolate_second_order_sw_fun" : 256,
        "extrapolate_velocity_second_order_true_fun" : 256,
        "extrapolate_second_order_edge_swb2_fun" : 256,
        "extrapolate_second_order_and_limit_by_vertex_fun" : 256,
        "extrapolate_second_order_and_limit_by_edge_fun" : 256, # FIXME

        "evaluate_segment_reflective_fun" : 256,
        "evaluate_segment_dirichlet_1_fun" : 256, # FIXME
        "evaluate_segment_dirichlet_2_fun" : 256, # FIXME

        "gravity_fun" : 448,

        "get_absolute_fun" : 256,

        "interpolate_fun" : 256,

        "manning_friction_flat_fun" : 256,
        "manning_friction_sloped_fun" : 256,

        "protect_sw_ext_fun" : 256,
        "protect_swb2_fun" : 256, # FIXME

        "saxpy_fun" : 256,

        "set_boundary_values_from_edges_fun" : 256,

        "update_centroids_fun" : 256,
        "update_fun" : 256,
        }
elif WorkingStation == "GTX480":
    kernel_block_configuration = {
        "balance_fun" : 256,

        "compute_fluxes_fun" : 128,

        "extrapolate_first_order_fun" : 1024,

        "extrapolate_second_order_sw_fun" : 128,
        "extrapolate_velocity_second_order_true_fun" : 128,
        "extrapolate_second_order_edge_swb2_fun" : 64,
        "extrapolate_second_order_and_limit_by_vertex_fun" : 256,
        "extrapolate_second_order_and_limit_by_edge_fun" : 256, # FIXME

        "evaluate_segment_reflective_fun" : 256,
        "evaluate_segment_dirichlet_1_fun" : 256, # FIXME
        "evaluate_segment_dirichlet_2_fun" : 256, # FIXME

        "gravity_fun" : 512,

        "get_absolute_fun" : 512,

        "interpolate_fun" : 1024,

        "manning_friction_flat_fun" : 192,
        "manning_friction_sloped_fun" : 128,

        "protect_sw_ext_fun" : 512,
        "protect_swb2_fun" : 64, # FIXME

        "saxpy_fun" : 512,

        "set_boundary_values_from_edges_fun" :512,

        "update_centroids_fun" : 512,
        "update_fun" : 512,
        }
else:
    kernel_block_configuration = {
        "balance_fun" : 32,

        "compute_fluxes_fun" : 64,

        "extrapolate_first_order_fun" : 32,

        "extrapolate_second_order_sw_fun" : 64,
        "extrapolate_velocity_second_order_true_fun" : 128,
        "extrapolate_second_order_edge_swb2_fun" : 64,
        "extrapolate_second_order_and_limit_by_vertex_fun" : 64,
        "extrapolate_second_order_and_limit_by_edge_fun" : 64, # FIXME

        "evaluate_segment_reflective_fun" : 64,
        "evaluate_segment_dirichlet_1_fun" : 64, # FIXME
        "evaluate_segment_dirichlet_2_fun" : 64, # FIXME

        "gravity_fun" :32,

        "get_absolute_fun" : 64,

        "interpolate_fun" : 32,

        "manning_friction_flat_fun" : 192,
        "manning_friction_sloped_fun" : 128,

        "protect_sw_ext_fun" : 64,
        "protect_swb2_fun" : 64, # FIXME

        "saxpy_fun" : 64,

        "set_boundary_values_from_edges_fun" : 32,

        "update_centroids_fun" : 64,
        "update_fun" : 64,
        }


balance_stream =                0
cf_central_stream =             1
extra_1_stream =                2
extra_2_sw_stream =             3
extra_velocity_2_stream =       4
extra_2_edge_swb2_stream =      5
extra_2_limit_vertex_stream =   6 
extra_2_limit_edge_stream =     7
evaluate_seg_ref_stream =       8
evaluate_seg_dir_1_stream =     9
evaluate_seg_dir_2_stream =     10
gravity_wb_stream =             11
interpolate_stream =            12
manning_flat_stream =           13
manning_sloped_stream =         14
protect_sw_stream =            15
protect_swb2_stream =           16
saxpy_stream =                  17
set_boundary_from_edge_stream = 18
update_centroid_stream =        19
update_stream =                 20
