#!/usr/bin/env python


import numpy
from pycuda import driver as drv
from anuga_cuda import generate_merimbula_domain
from anuga_cuda import generate_cairns_domain

using_tsunami_domain = True

if using_tsunami_domain:
    domain1 = generate_cairns_domain(False)
    domain2 = generate_cairns_domain(True)
else:
    domain1 = generate_merimbula_domain()
    domain2 = generate_merimbula_domain(gpu=True)

yieldstep = 50
finaltime = 500

domain1.evolve(yieldstep = yieldstep, finaltime = finaltime)
domain1.evolve(yieldstep = yieldstep, finaltime = finaltime)
domain1.evolve(yieldstep = yieldstep, finaltime = finaltime)
    

domain2.evolve(yieldstep = yieldstep, finaltime = finaltime)
domain2.evolve(yieldstep = yieldstep, finaltime = finaltime)
domain2.evolve(yieldstep = yieldstep, finaltime = finaltime)

N = domain2.number_of_elements
W1 = 32
W2 = 1
W3 = 1

#domain1.yieldstep = domain1.get_time() + yieldstep
#domain1.evolve_one_rk2_step(yieldstep = yieldstep, finaltime=finaltime)
domain1.backup_conserved_quantities()

#domain2.yieldstep = domain2.get_time() + yieldstep
domain2.backup_conserved_quantities()
#domain2.compute_fluxes()
#domain2.compute_forcing_terms()
#domain2.update_timestep(yieldstep = yieldstep, finaltime = finaltime)
#domain2.update_conserved_quantities()
#domain2.set_time(domain2.get_time() + domain2.timestep)
#domain2.update_boundary()
#domain2.compute_fluxes()
#domain2.compute_forcing_terms()
#domain2.update_conserved_quantities()

a = 0.5
b = 0.5

for name in domain2.conserved_quantities:
    Q1 = domain1.quantities[name]
    Q2 = domain2.quantities[name]
    c1 = Q1.centroid_values
    c2 = Q2.centroid_values
    if not  numpy.allclose(c1, c2 ):
        print c1
        print c2 
    Q1.saxpy_centroid_values(a, b)
    domain2.saxpy_centroid_values_func(
            numpy.int32( N ),
            numpy.float64( a ),
            numpy.float64( b ),
            drv.InOut( Q2.centroid_values ),
            drv.In( Q2.centroid_backup_values ),
            block = (W1, W2, W3),
            grid=(( N +W1*W2*W3-1)/(W1*W2*W3),1)
            )




    res = numpy.allclose(c1, c2) 

    print "\n\n centroid_values all closed? ", res

    if not res:
        print c1
        print c2
        cnt_c = 0
        for i in range(N):
            if (c1[i] != c2[i]).all():
                cnt_c += 1
                if cnt_c <= 5:
                    print i, c1[i], c2[i]
        print " Number of diff: is %d " % cnt_c
