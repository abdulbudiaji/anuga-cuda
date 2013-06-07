#!/usr/bin/env python


import numpy
from pycuda import driver as drv
from anuga_cuda import generate_merimbula_domain, generate_cairns_domain, get_kernel_function_info

using_tsunami_domain = False

if using_tsunami_domain:
    domain1 = generate_cairns_domain(False)
    domain2 = generate_cairns_domain(True)
else:
    domain1 = generate_merimbula_domain()
    domain2 = generate_merimbula_domain(gpu=True)

domain2.equip_kernel_functions()

domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
if using_tsunami_domain:
    domain1.distribute_to_vertices_and_edges()
else:
    from swb2_domain_ext import \
            extrapolate_second_order_edge_sw as extrapol2_ext

    # Shortcuts
    Stage = domain1.quantities['stage']
    Xmom = domain1.quantities['xmomentum']
    Ymom = domain1.quantities['ymomentum']
    Elevation = domain1.quantities['elevation']
    
    extrapol2_ext(domain1,
            domain1.surrogate_neighbours,
            domain1.number_of_boundaries,
            domain1.centroid_coordinates,
            Stage.centroid_values,
            Xmom.centroid_values,
            Ymom.centroid_values,
            Elevation.centroid_values,
            domain1.edge_coordinates,
            Stage.edge_values,
            Xmom.edge_values,
            Ymom.edge_values,
            Elevation.edge_values,
            Stage.vertex_values,
            Xmom.vertex_values,
            Ymom.vertex_values,
            Elevation.vertex_values,
            int(domain1.optimise_dry_cells),
            int(domain1.extrapolate_velocity_second_order))


N = domain2.number_of_elements
domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)
if using_tsunami_domain:
    domain2.protect_against_infinitesimal_and_negative_heights()

import sys
W1 = 0
for i in range( len(sys.argv)):
    if sys.argv[i] == "-b":
        W1 = int(sys.argv[i+1])

if not W1:
    W1 = domain2.extrapolate_second_order_edge_swb2_func.max_threads_per_block
print W1
W2 = 1
W3 = 1

get_kernel_function_info(
    domain2.extrapolate_second_order_edge_swb2_func, W1,W2, W3)

domain2.stage_centroid_store = numpy.zeros(N, dtype=numpy.float64)
domain2.xmomentum_centroid_store = numpy.zeros(N, dtype=numpy.float64)
domain2.ymomentum_centroid_store = numpy.zeros(N, dtype=numpy.float64)
domain2.mim_elevation_edgevalue = numpy.zeros(N, dtype=numpy.float64)
domain2.max_elevation_edgevalue = numpy.zeros(N, dtype=numpy.float64)
domain2.count_wet_neighbours = numpy.zeros(N, dtype=numpy.int32)

domain2.extrapolate_second_order_edge_swb2_func(
        numpy.int32(N),
        numpy.int32(domain2.optimise_dry_cells),
        numpy.int32(domain2.extrapolate_velocity_second_order),
        numpy.float64(domain2.epsilon),
        numpy.float64(domain2.minimum_allowed_height),
        numpy.float64(domain2.beta_w),
        numpy.float64(domain2.beta_w_dry),
        numpy.float64(domain2.beta_uh),
        numpy.float64(domain2.beta_uh_dry),
        numpy.float64(domain2.beta_vh),
        numpy.float64(domain2.beta_vh_dry),

        drv.In( domain2.surrogate_neighbours),
        drv.In( domain2.number_of_boundaries),
        drv.In( domain2.centroid_coordinates),
        
        drv.In( domain2.quantities['stage'].centroid_values),
        drv.In( domain2.quantities['elevation'].centroid_values),
        drv.InOut( domain2.quantities['xmomentum'].centroid_values),
        drv.InOut( domain2.quantities['ymomentum'].centroid_values),

        drv.In( domain2.edge_coordinates),

        drv.InOut( domain2.quantities['stage'].edge_values),
        drv.In( domain2.quantities['elevation'].edge_values),
        drv.InOut( domain2.quantities['xmomentum'].edge_values),
        drv.InOut( domain2.quantities['ymomentum'].edge_values),

        drv.InOut( domain2.quantities['stage'].vertex_values),
        drv.In( domain2.quantities['elevation'].vertex_values),
        drv.InOut( domain2.quantities['xmomentum'].vertex_values),
        drv.InOut( domain2.quantities['ymomentum'].vertex_values),
        
        drv.In( domain2.stage_centroid_store),
        drv.In( domain2.xmomentum_centroid_store),
        drv.In( domain2.ymomentum_centroid_store),

        drv.InOut( domain2.mim_elevation_edgevalue),
        drv.InOut( domain2.max_elevation_edgevalue),
        drv.InOut( domain2.count_wet_neighbours),

        block = (W1, W2, W3),
        grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
    )


xc1 = domain1.quantities['xmomentum'].centroid_values
yc1 = domain1.quantities['ymomentum'].centroid_values
se1 = domain1.quantities['stage'].edge_values
xe1 = domain1.quantities['xmomentum'].edge_values
ye1 = domain1.quantities['ymomentum'].edge_values
sv1 = domain1.quantities['stage'].vertex_values
xv1 = domain1.quantities['xmomentum'].vertex_values
yv1 = domain1.quantities['ymomentum'].vertex_values

xc2 = domain2.quantities['xmomentum'].centroid_values
yc2 = domain2.quantities['ymomentum'].centroid_values
se2 = domain2.quantities['stage'].edge_values
xe2 = domain2.quantities['xmomentum'].edge_values
ye2 = domain2.quantities['ymomentum'].edge_values
sv2 = domain2.quantities['stage'].vertex_values
xv2 = domain2.quantities['xmomentum'].vertex_values
yv2 = domain2.quantities['ymomentum'].vertex_values

res = []
res.append( numpy.allclose(xc1, xc2) )
res.append( numpy.allclose(yc1, yc2) )
res.append( numpy.allclose(se1, se2) )
res.append( numpy.allclose(xe1, xe2) )
res.append( numpy.allclose(ye1, ye2) )
res.append( numpy.allclose(sv1, sv2) )
res.append( numpy.allclose(xv1, xv2) )
res.append( numpy.allclose(yv1, yv2) )

print "\n\nxmom_centroid_values all closed? " ,res[0]
print "ymom_centroid_values all closed? " ,res[1]
print "stage_edge_values all closed?  " ,res[2]
print "xmom_edge_values all closed?  " ,res[3]
print "ymom_edge_values all closed?  " ,res[4]
print "stage_vertex_values all closed?  " ,res[5]
print "xmom_vertex_values all closed?  " ,res[6]
print "ymom_vertex_values all closed?  " ,res[7]

if not res.count(True) == 8:
    cnt_xc = 0
    cnt_yc = 0
    cnt_se = 0
    cnt_xe = 0
    cnt_ye = 0
    cnt_sv = 0
    cnt_xv = 0
    cnt_yv = 0
    res = []
    for i in range(N):
        if (se1[i] != se2[i]).any():
            cnt_se += 1
            res.append(i)
            if cnt_se :
                print i, se1[i], se2[i]
    print " Number of diff: is %d " % cnt_se
    print res
