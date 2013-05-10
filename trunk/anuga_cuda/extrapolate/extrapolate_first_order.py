#!/usr/bin/env python

import numpy
from pycuda import driver as drv
from anuga_cuda import *

testing_rearranged = False

domain1 = generate_merimbula_domain( gpu=False )
domain2 = generate_merimbula_domain( gpu=True )
if testing_rearranged:
    domain2 = rearrange_domain(domain2)
    sort_domain(domain1)
domain2.equip_kernel_functions()

domain1.protect_against_infinitesimal_and_negative_heights()
domain2.protect_against_infinitesimal_and_negative_heights()

N = domain1.number_of_elements
import sys
W1 = 0
for i in range( len(sys.argv)):
    if sys.argv[i] == "-b":
        W1 = int(sys.argv[i+1])

if not W1:
    W1 = domain2.extrapolate_first_order_func.max_threads_per_block
print W1
W2 = 1
W3 = 1

get_kernel_function_info(domain2.extrapolate_first_order_func,W1,W2, W3)


for name in domain1.conserved_quantities:
    Q1 = domain1.quantities[name]
    Q1.extrapolate_first_order()

    Q2 = domain2.quantities[name]
    domain2.extrapolate_first_order_func(
        numpy.int32(N),
        drv.In( Q2.centroid_values ),
        drv.InOut( Q2.edge_values ),
        drv.InOut( Q2.vertex_values ),
        block = (W1, W2, W3),
        grid = ( (N +W1*W2*W3 -1)/(W1*W2*W3), 1)
        )

    cnt_e = 0
    cnt_v = 0

    print "*********** in %s\n" % name
    e1 = Q1.edge_values
    v1 = Q1.vertex_values
    e2 = Q2.edge_values
    v2 = Q2.vertex_values
    for i in range(N):
        if (Q1.edge_values[i] != Q2.edge_values[i]).all():
        #if e1[i][0] != e2[i/3][i%3] or \
        #        e1[i][1] != e2[(i+N)/3][(i+N)%3] or \
        #        e1[i][2] != e2[(i+2*N)/3][(i+2*N)%3]:
            cnt_e += 1
            if cnt_e  < 5:
                print i, Q1.edge_values[i], Q2.edge_values[i]
        if (Q1.vertex_values[i] != Q2.vertex_values[i]).all():
        #if v1[i][0] != v2[i/3][i%3] or \
        #        v1[i][1] != v2[(i+N)/3][(i+N)%3] or \
        #        v1[i][2] != v2[(i+2*N)/3][(i+2*N)%3]:
            cnt_v += 1
            if cnt_v  < 5:
                print i, Q1.vertex_values[i], Q2.vertex_values[i]

    print "___________ %d, %d\n" % ( cnt_e, cnt_v)
