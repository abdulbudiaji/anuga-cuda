#!/usr/bin/env python

import numpy 
from utilities import cpy_back_and_cmp



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


def approx_equal(a, b, approx=True):
    """Check result with tolerance.
    
    return abs(a-b) <= abs(a)*tolerance
    """

    if approx:
        if abs(a-b) > abs(a)*pow(10,-6):
            return False
        else: 
            return True
    else:
        if a != b:
            return False
        else:
            return True





def check_result(a, b):
    """Check result correctness. """

    global global_cnt
    global name_list
    if not numpy.allclose(a, b):
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
    """Check all the result"""

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



def test_distribute_to_vertexs_and_edges(domain, IO = 'Output'):
    """Pair-testing check point for distribute_to_vertices_and_edges """

    # For time-dependence issues, used to synchronize the kernel
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

    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> distribute_to_vertices_and_edges ", res)




def test_evolve_one_euler_step(domain):
    """Pari-testing check point for evolve_one_euler_step"""

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
    

    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> evolve_one_euler_step ", res)




def test_update_ghosts(domain):
    """Pari-testing check point for update_ghosts"""

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


    if res.count(True) + res.count(-1) != res.__len__():
        print " --> update_ghosts ",res




def test_update_extrema(domain):
    """Pari-testing check point for update_extrema"""

    # For time-dependence issues
    #ctx.synchronize()

    gpu = domain.using_gpu
    rg = domain.rearranged_domain
    sc = domain.cotesting_domain

    res = []
    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> update_extrema ", res)
    


def test_update_timestep(domain):
    """Pari-testing check point for update_timestep"""

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


    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> update_timestep ",res)



def test_update_conserved_quantities(domain, output =True):
    """Pari-testing check point for update_conserved_quantities"""

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
        
        if res.count(True) + res.count(-1) != res.__len__():
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
    """Pari-testing check point for manning_friction_implicit"""

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

    if res.count(True) + res.count(-1) != res.__len__():
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
    """Pari-testing check point for update_boundary"""

    # For time-dependence issues
    #ctx.synchronize()

    from anuga.shallow_water.boundaries import Reflective_boundary
    from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
            import Transmissive_boundary, Dirichlet_boundary, \
                    Compute_fluxes_boundary, Time_boundary, Time_boundary,\
                    Time_boundary


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

            if ipt.count(True) + ipt.count(0) != ipt.__len__():
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
    """Pari-testing check point for update_other_quantities"""

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

    if res.count(True) + res.count(-1) != res.__len__():
       raise Exception("Error in update_other_quantities", res)



def test_compute_fluxes(domain):
    """Pari-testing check point for compute_fluxes"""

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
    

    if res.count(True) + res.count(-1) != res.__len__():
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
    """Pari-testing check point for compute_forcing_terms"""

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

    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> compute_forcing_terms ", res)



def test_protect_against_infinitesimal_and_negative_heights(domain):
    """Pari-testing check point for 
        protect_against_infinitesimal_and_negative_heights
    """

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


    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> protect_against_infinitesimal_and_negative_heights ",res)



def test_extrapolate_second_order_sw(domain):
    """Pari-testing check point for extrapolate_second_order_sw"""

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

    if res.count(True) + res.count(-1) != res.__len__():
        print res

        if False:
            for i in range(domain.number_of_elements):
                if (xm1.vertex_values[i] != xm2.vertex_values[i]).all():
                    if domain.number_of_boundaries[i] == 1:
                        print i, xm1.vertex_values[i], xm2.vertex_values[i]
                    cnt += 1
            print cnt

        raise Exception( " --> extrapolate_second_order_sw ", res)



def test_balance_deep_and_shallow(domain):
    """Pari-testing check point for balance_deep_and_shallow"""

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


    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> balance_deep_and_shallow ", res)
        


def test_interpolate_from_vertices_to_edges(domain):
    """Pari-testing check point for interpolate_from_vertices_to_edges"""

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
    res.append( cpy_back_and_cmp( s1, s2, 'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( s1, s2, 'edge_values', gpu, rg))

    res.append( cpy_back_and_cmp( xm1, xm2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( xm1, xm2,'edge_values', gpu, rg))

    res.append( cpy_back_and_cmp( ym1, ym2,'vertex_values', gpu, rg))
    res.append( cpy_back_and_cmp( ym1, ym2,'edge_values', gpu, rg))


    if res.count(True) + res.count(-1) != res.__len__():
        for i in range(domain.number_of_elements):
            if not numpy.allclose(s1.edge_values[i], s2.edge_values[i]):
                print i, s1.edge_values[i], s2.edge_values[i]
        raise Exception( " --> interpolate_from_vertices_to_edges ", res)



def test_extrapolate_second_order_and_limit_by_vertex(domain):
    """Pari-testing check point for 
        extrapolate_second_order_and_limit_by_vertex
    """

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


    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> extrapolate_second_order_and_limit_by_vertex ", res)

    

def test_extrapolate_first_order(domain):
    """Pari-testing check point for extrapolate_first_order"""

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


    if res.count(True) + res.count(-1) != res.__len__():
        raise Exception( " --> extrapolate_first_order ", res)



def test_update_centroids_of_velocities_and_height(domain):
    """Pari-testing check point for 
        update_centroids_of_velocities_and_height
    """

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

    if res.count(True) + res.count(-1) != res.__len__():
        if not res[12]:
            for i in range(domain.number_of_elements):
                if not numpy.allclose(
                        xv1.centroid_values[i],
                        xv2.centroid_values[i]) :
                    print "  xvelocity centroid %d %lf %lf" % \
                        (i, xv1.centroid_values[i], xv2.centroid_values[i])
        raise Exception("Error in update_centroids_of_velocities_and_height", res)


