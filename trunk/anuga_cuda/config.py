testing_domain_path="/home/u5044856/anuga-cuda/anuga_cuda/testing_domains/"

merimbula_dir = testing_domain_path+"merimbula/"

cairns_dir= testing_domain_path+"cairns/"


kernel_path = {
    "compute_fluxes_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/compute_fluxes/",
    "gravity_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/gravity/",
    "extrapolate_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/extrapolate/",
    "protect_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/protect/",
    "balance_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/balance/",
    "interpolate_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/interpolate_from_vertices_to_edges/",
    "evaluate_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/evaluate_segment/",
    "get_absolute_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/get_absolute/",
    "manning_friction_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/manning_friction/",
    "saxpy_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/saxpy_centroid_values/",
    "set_boundary_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/set_boundary/",
    "update_centroids_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/update_centroids_of_velocities_and_height/",
    "update_dir" : \
        "/home/u5044856/anuga-cuda/anuga_cuda/update/"
    }


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

    "gravity_fun" : 32,

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
