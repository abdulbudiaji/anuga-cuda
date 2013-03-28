#!/usr/bin/env python


import numpy
from pycuda import driver as drv
from anuga_cuda import *


using_tsunami_domain = False
using_rearranged_domain = True

if using_tsunami_domain:
    domain1 = generate_cairns_domain(False)
    domain2 = generate_cairns_domain(True)
else:
    domain1 = generate_merimbula_domain()
    domain2 = generate_merimbula_domain(gpu=True)

if using_rearranged_domain:
    domain2 = rearrange_domain(domain2)
    sort_domain(domain1)


domain2.equip_kernel_functions()
import sys
W1 = 0
for i in range( len(sys.argv)):
    if sys.argv[i] == "-b":
        W1 = int(sys.argv[i+1])

if not W1:
    W1 = domain2.update_func.max_threads_per_block
W2 = 1
W3 = 1

get_kernel_function_info(domain2.update_func, 
        W1,W2, W3)


domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
    

domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)

N = domain2.number_of_elements


domain1.update_conserved_quantities()

for name in domain2.conserved_quantities:
    Q2 = domain2.quantities[name]

    domain2.update_func(
        numpy.int32(N),
        numpy.float64(domain2.timestep),

        drv.InOut( Q2.centroid_values ),
        drv.In( Q2.explicit_update ),

        drv.InOut( Q2.semi_implicit_update ),

        block = (W1, W2, W3),
        grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
        )
        
    Q1 = domain1.quantities[name]
    c1 = Q1.centroid_values
    c2 = Q2.centroid_values
    s1 = Q1.semi_implicit_update
    s2 = Q2.semi_implicit_update


    res = []
    res.append( numpy.allclose(c1, c2) )
    res.append( numpy.allclose(s1, s2) )

    print "%s_centroid_values all closed? %s" % \
            (name,'True' if res[0] else 'False')
    print "%s_semi_implicit_update all closed? %s" % \
            (name, 'True' if res[1] else 'False')

    if not res.count(True) == 2:
        cnt_c = 0
        cnt_s = 0
        res = []
        for i in range(N):
            if (s1[i] != s2[i]).any():
                cnt_s += 1
                res.append(i)
                if cnt_s :
                    print i, s1[i], s2[i]
        print " Number of diff: is %d " % cnt_s
        print res
