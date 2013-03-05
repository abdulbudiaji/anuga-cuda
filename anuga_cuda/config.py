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
    "compute_fluxes_fun" : 64,

    "gravity_fun" : 32,

    "extrapolate_first_order_fun" : 32,

    "extrapolate_second_order_sw_fun" : 64,
    "extrapolate_velocity_second_order_true_fun" : 128,
    "extrapolate_second_order_edge_swb2_fun" : 64,
    "extrapolate_second_order_and_limit_by_vertex_fun" : 64,
    "extrapolate_second_order_and_limit_by_edge_fun" : 64, # FIXME

    "protect_sw_ext_fun" : 64,
    "protect_swb2_fun" : 64, # FIXME

    "balance_fun" : 32,

    "interpolate_fun" : 32,

    "evaluate_segment_reflective_fun" : 64,
    "evaluate_segment_dirichlet_1_fun" : 64, # FIXME
    "evaluate_segment_dirichlet_2_fun" : 64, # FIXME

    "get_absolute_fun" : 64,

    "manning_friction_flat_fun" : 192,
    "manning_friction_sloped_fun" : 128,

    "saxpy_fun" : 64,

    "set_boundary_values_from_edges_fun" : 32,

    "update_centroids_fun" : 64,
    "update_fun" : 64,
    }

