#!/usr/bin/env python


from utilities.utility import get_sourceModule, get_page_locked_array, \
        get_kernel_function_info, get_device_array, asy_cpy, cpy_back, \
        number_domain_method, cpy_back_and_cmp

from utilities.check_result import approx_equal, check_result, check_all, \
        test_distribute_to_vertexs_and_edges, test_evolve_one_euler_step,\
        test_update_ghosts, test_update_extrema, test_update_timestep, \
        test_update_conserved_quantities, test_manning_friction_implicit, \
        test_update_boundary, test_update_other_quantities, \
        test_update_centroids_of_velocities_and_height, \
        test_compute_fluxes, test_compute_forcing_terms, \
        test_protect_against_infinitesimal_and_negative_heights,\
        test_extrapolate_second_order_sw, test_balance_deep_and_shallow,\
        test_interpolate_from_vertices_to_edges, \
        test_extrapolate_second_order_and_limit_by_vertex
        

from utilities.sort_domain import swap_domain, sort_domain, \
        sort_domain_check, rearrange_domain, check_rearranged_array, \
        rearrange_domain_check
