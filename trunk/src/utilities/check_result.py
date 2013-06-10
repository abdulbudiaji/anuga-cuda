#!/usr/bin/env python


import numpy as num



global_cnt = 0
name_list = ['max_speed', 
                'stage_edge', 'elevation_edge', 'xmom_edge', 'ymom_edge', 
                'height_edge', 'xvelocity_edge', 'yvelocity_edge', 
                
                'stage_centroid', 'elevation_centroid', 'xmom_centroid', 'ymom_centroid', 
                'height_centroid', 'xvelocity_centroid', 'yvelocity_centroid', 
                'friction_centroid',
                
                'stage_vertex', 'elevation_vertex', 'xmom_vertex', 'ymom_vertex', 
                'height_vertex', 'xvelocity_vertex', 'yvelocity_vertex', 
                
                'stage_explicit_update', 'xmom_explicit_update', 'ymom_explicit_update', 
                'stage_semi_implicit_update', 'xmom_semi_implicit_update', 
                'ymom_semi_implicit_update', 

                'stage_x_gradient', 'elevation_x_gradient', 'xmom_x_gradient', 
                'ymom_x_gradient', 'height_x_gradient', 'xvelocity_x_gradient', 
                'yvelocity_x_gradient', 
                
                'stage_y_gradient', 'elevation_y_gradient', 'xmom_y_gradient', 
                'ymom_y_gradient', 'height_y_gradient', 'xvelocity_y_gradient', 
                'yvelocity_y_gradient', 

                'stage_boundary', 'elevation_boundary', 'xmom_boundary', 'ymom_boundary', 
                'height_boundary', 'xvelocity_boundary', 'yvelocity_boundary' 
            ]



def check_result(a, b):
    global global_cnt
    global name_list
    if not num.allclose(a, b):
        print global_cnt, name_list[global_cnt]
        cnt = 0
        if a.shape.__len__()  == 1:
            for i in range(a.shape[0]):
                if ( a[i] != b[i]):
                    cnt += 1
                    if cnt <= 10:
                        print i, a[i], b[i] 
        else:
            for i in range(a.shape[0]):
                if ( a[i] != b[i]).any():
                    cnt += 1
                    if cnt <= 10:
                        print i, a[i], b[i] 
        print cnt
    global_cnt += 1



def check_all(domain, test_domain):
    global global_cnt
    global_cnt = 0
    print domain.flux_timestep, test_domain.flux_timestep
    
    
    check_result(domain.max_speed, test_domain.max_speed)
    
    
    check_result(domain.quantities['stage'].edge_values, 
            test_domain.quantities['stage'].edge_values)
    check_result(domain.quantities['elevation'].edge_values, 
            test_domain.quantities['elevation'].edge_values)
    check_result(domain.quantities['xmomentum'].edge_values, 
            test_domain.quantities['xmomentum'].edge_values)
    check_result(domain.quantities['ymomentum'].edge_values, 
            test_domain.quantities['ymomentum'].edge_values)
    check_result(domain.quantities['height'].edge_values, 
            test_domain.quantities['height'].edge_values)
    check_result(domain.quantities['xvelocity'].edge_values, 
            test_domain.quantities['xvelocity'].edge_values)
    check_result(domain.quantities['yvelocity'].edge_values, 
            test_domain.quantities['yvelocity'].edge_values)
    
    
    check_result(domain.quantities['stage'].centroid_values, 
            test_domain.quantities['stage'].centroid_values)
    check_result(domain.quantities['elevation'].centroid_values, 
            test_domain.quantities['elevation'].centroid_values)
    check_result(domain.quantities['xmomentum'].centroid_values, 
            test_domain.quantities['xmomentum'].centroid_values)
    check_result(domain.quantities['ymomentum'].centroid_values, 
            test_domain.quantities['ymomentum'].centroid_values)
    check_result(domain.quantities['height'].centroid_values, 
            test_domain.quantities['height'].centroid_values)
    check_result(domain.quantities['xvelocity'].centroid_values, 
            test_domain.quantities['xvelocity'].centroid_values)
    check_result(domain.quantities['yvelocity'].centroid_values, 
            test_domain.quantities['yvelocity'].centroid_values)
    check_result(domain.quantities['friction'].centroid_values, 
            test_domain.quantities['friction'].centroid_values)
    
    
    check_result(domain.quantities['stage'].vertex_values, 
            test_domain.quantities['stage'].vertex_values)
    check_result(domain.quantities['elevation'].vertex_values, 
            test_domain.quantities['elevation'].vertex_values)
    check_result(domain.quantities['xmomentum'].vertex_values, 
            test_domain.quantities['xmomentum'].vertex_values)
    check_result(domain.quantities['ymomentum'].vertex_values, 
            test_domain.quantities['ymomentum'].vertex_values)
    check_result(domain.quantities['height'].vertex_values, 
            test_domain.quantities['height'].vertex_values)
    check_result(domain.quantities['xvelocity'].vertex_values, 
            test_domain.quantities['xvelocity'].vertex_values)
    check_result(domain.quantities['yvelocity'].vertex_values, 
            test_domain.quantities['yvelocity'].vertex_values)
    
    
    check_result(domain.quantities['stage'].explicit_update, 
            test_domain.quantities['stage'].explicit_update)
    check_result(domain.quantities['xmomentum'].explicit_update, 
            test_domain.quantities['xmomentum'].explicit_update)
    check_result(domain.quantities['ymomentum'].explicit_update, 
            test_domain.quantities['ymomentum'].explicit_update)
    
        
    check_result(domain.quantities['stage'].semi_implicit_update, 
            test_domain.quantities['stage'].semi_implicit_update)
    check_result(domain.quantities['xmomentum'].semi_implicit_update, 
            test_domain.quantities['xmomentum'].semi_implicit_update)
    check_result(domain.quantities['ymomentum'].semi_implicit_update, 
            test_domain.quantities['ymomentum'].semi_implicit_update)
    
    
    check_result(domain.quantities['stage'].x_gradient, 
            test_domain.quantities['stage'].x_gradient)
    check_result(domain.quantities['elevation'].x_gradient, 
            test_domain.quantities['elevation'].x_gradient)
    check_result(domain.quantities['xmomentum'].x_gradient, 
            test_domain.quantities['xmomentum'].x_gradient)
    check_result(domain.quantities['ymomentum'].x_gradient, 
            test_domain.quantities['ymomentum'].x_gradient)
    check_result(domain.quantities['height'].x_gradient, 
            test_domain.quantities['height'].x_gradient)
    check_result(domain.quantities['xvelocity'].x_gradient, 
            test_domain.quantities['xvelocity'].x_gradient)
    check_result(domain.quantities['yvelocity'].x_gradient, 
            test_domain.quantities['yvelocity'].x_gradient)
    
    
    check_result(domain.quantities['stage'].y_gradient, 
            test_domain.quantities['stage'].y_gradient)
    check_result(domain.quantities['elevation'].y_gradient, 
            test_domain.quantities['elevation'].y_gradient)
    check_result(domain.quantities['xmomentum'].y_gradient, 
            test_domain.quantities['xmomentum'].y_gradient)
    check_result(domain.quantities['ymomentum'].y_gradient, 
            test_domain.quantities['ymomentum'].y_gradient)
    check_result(domain.quantities['height'].y_gradient, 
            test_domain.quantities['height'].y_gradient)
    check_result(domain.quantities['xvelocity'].y_gradient, 
            test_domain.quantities['xvelocity'].y_gradient)
    check_result(domain.quantities['yvelocity'].y_gradient, 
            test_domain.quantities['yvelocity'].y_gradient)
    
    
    check_result(domain.quantities['stage'].boundary_values, 
            test_domain.quantities['stage'].boundary_values)
    check_result(domain.quantities['elevation'].boundary_values, 
            test_domain.quantities['elevation'].boundary_values)
    check_result(domain.quantities['xmomentum'].boundary_values, 
            test_domain.quantities['xmomentum'].boundary_values)
    check_result(domain.quantities['ymomentum'].boundary_values, 
            test_domain.quantities['ymomentum'].boundary_values)
    check_result(domain.quantities['height'].boundary_values, 
            test_domain.quantities['height'].boundary_values)
    check_result(domain.quantities['xvelocity'].boundary_values, 
            test_domain.quantities['xvelocity'].boundary_values)
    check_result(domain.quantities['yvelocity'].boundary_values, 
            test_domain.quantities['yvelocity'].boundary_values)



def number_domain_method(domain):
    if domain.timestepping_method == 'euler':
        timestepping_method = 1
    elif domain.timestepping_method == 'rk2':
        timestepping_method = 2
    elif domain.timestepping_method == 'rk3':
        timestepping_method = 3
    else:
        timestepping_method = 4
    print " The timestepping_method is '%s' %d" % (domain.timestepping_method, timestepping_method)
    
    
    if domain.flow_algorithm == 'tsunami':
        flow_algorithm = 1
    elif domain.flow_algorithm == 'yusuke':
        flow_algorithm = 2
    else:
        flow_algorithm = 3
    print " The flow_algorithm us '%s' %d" % (domain.flow_algorithm, flow_algorithm)
    
    
    if domain.compute_fluxes_method == 'original':
        compute_fluxes_method = 0
    elif domain.compute_fluxes_method == 'wb_1':
        compute_fluxes_method = 1
    elif domain.compute_fluxes_method == 'wb_2':
        compute_fluxes_method = 2
    elif domain.compute_fluxes_method == 'wb_3':
        compute_fluxes_method = 3
    elif domain.compute_fluxes_method == 'tsunami':
        compute_fluxes_method = 4
    else:
        compute_fluxes_method = 5
    print " The compute_fluxes_method is '%s' %d" % (domain.compute_fluxes_method, compute_fluxes_method)


    return (compute_fluxes_method, flow_algorithm, timestepping_method)
