#!/usr/bin/env python


import numpy
from pycuda import driver as drv
from anuga_cuda import generate_merimbula_domain
from anuga_cuda import generate_cairns_domain

#domain1 = generate_merimbula_domain()
#domain2 = generate_merimbula_domain(gpu=True)
domain1 = generate_cairns_domain(False)
domain2 = generate_cairns_domain(True)

domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.protect_against_infinitesimal_and_negative_heights()

N = domain2.number_of_elements
domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)

W1 = 32
W2 = 1
W3 = 1
domain2.protect_swb2_func(
        numpy.int32(N),
        numpy.float64(domain2.minimum_allowed_height),
        numpy.float64(domain2.maximum_allowed_speed),
        numpy.float64(domain2.epsilon),
        drv.InOut(domain2.quantities['stage'].centroid_values),
        drv.InOut(domain2.quantities['stage'].vertex_values),
        drv.In(domain2.quantities['elevation'].centroid_values),
        drv.In(domain2.quantities['elevation'].vertex_values),
        drv.InOut(domain2.quantities['xmomentum'].centroid_values),
        drv.InOut(domain2.quantities['ymomentum'].centroid_values),
        drv.In(domain2.areas),
        block = (W1, W2, W3),
        grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
        )


cnt_sc = 0
cnt_sv = 0
cnt_ec = 0
cnt_ev = 0
cnt_xc = 0
cnt_yc = 0

sc1 = domain1.quantities['stage'].centroid_values
sv1 = domain1.quantities['stage'].vertex_values
xc1 = domain1.quantities['xmomentum'].centroid_values
yc1 = domain1.quantities['ymomentum'].centroid_values

sc2 = domain2.quantities['stage'].centroid_values
sv2 = domain2.quantities['stage'].vertex_values
xc2 = domain2.quantities['xmomentum'].centroid_values
yc2 = domain2.quantities['ymomentum'].centroid_values

print "stage_centroid_values all closed?" ,numpy.allclose(sc1, sc2)
print "stage_vertex_values all closed?  " ,numpy.allclose(sv1, sv2)
print "xmom_centroid_values all closed? " ,numpy.allclose(xc1, xc2)
print "ymom_centroid_values all closed? " ,numpy.allclose(yc1, yc2)

##for i in range(N):
    
