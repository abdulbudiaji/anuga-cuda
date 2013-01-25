#!/usr/bin/env python


import numpy
from pycuda import driver as drv
from anuga_cuda import generate_merimbula_domain
from anuga_cuda import generate_cairns_domain

using_tsunami_domain = False
testing_sloped=False

if using_tsunami_domain:
    domain1 = generate_cairns_domain(False)
    domain2 = generate_cairns_domain(True)
else:
    domain1 = generate_merimbula_domain()
    domain2 = generate_merimbula_domain(gpu=True)


domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
    

domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)

N = domain2.number_of_elements
W1 = 32
W2 = 1
W3 = 1


if testing_sloped:
    domain1.use_sloped_mannings = True
        
    domain2.manning_friction_sloped_func(
        numpy.int32(N),
        numpy.float64(domain2.g),
        numpy.float64(domain2.minimum_allowed_height),

        drv.In( domain2.vertex_coordinates ),
        drv.In( domain2.quantities['stage'].centroid_values ),
        drv.In( domain2.quantities['xmomentum'].centroid_values ),
        drv.In( domain2.quantities['ymomentum'].centroid_values ),
        drv.In( domain2.quantities['elevation'].vertex_values ),
        drv.In( domain2.quantities['friction'].centroid_values ),

        drv.InOut( domain2.quantities['xmomentum'].semi_implicit_update ),
        drv.InOut( domain2.quantities['ymomentum'].semi_implicit_update ),

        block = (W1, W2, W3),
        grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
        )
else:
    domain1.use_sloped_mannings = False
    domain2.manning_friction_flat_func(
            numpy.int32(N),
            numpy.float64(domain2.g),
            numpy.float64(domain2.minimum_allowed_height),

            drv.In( domain2.quantities['stage'].centroid_values ),
            drv.In( domain2.quantities['xmomentum'].centroid_values ),
            drv.In( domain2.quantities['ymomentum'].centroid_values ),
            drv.In( domain2.quantities['elevation'].centroid_values ),
            drv.In( domain2.quantities['friction'].centroid_values ),

            drv.In( domain2.quantities['xmomentum'].semi_implicit_update ),
            drv.In( domain2.quantities['ymomentum'].semi_implicit_update ),

            block = (W1, W2, W3),
            grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
            )
        
from anuga.shallow_water.shallow_water_domain import \
    manning_friction_implicit

manning_friction_implicit(domain1)



xs1 = domain1.quantities['xmomentum'].semi_implicit_update
ys1 = domain1.quantities['ymomentum'].semi_implicit_update

xs2 = domain2.quantities['xmomentum'].semi_implicit_update
ys2 = domain2.quantities['ymomentum'].semi_implicit_update


res = []
res.append( numpy.allclose(xs1, xs2) )
res.append( numpy.allclose(ys1, ys2) )

print "\n\nxmom_semi_implicit_update all closed? " ,res[0]
print "ymom_semi_implicit_update all closed? " ,res[1]

if not res.count(True) == 2:
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
