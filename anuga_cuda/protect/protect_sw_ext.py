#!/usr/bin/env python


import numpy
from pycuda import driver as drv
from anuga_cuda import *

using_rearranged_domain=True

domain1 = generate_merimbula_domain()
domain2 = generate_merimbula_domain(gpu=True)

if using_rearranged_domain:
    domain2 = rearrange_domain(domain2)
    sort_domain(domain1)


domain1.protect_against_infinitesimal_and_negative_heights()


domain2.equip_kernel_functions()
N = domain2.number_of_elements
import sys
W1 = 0
for i in range( len(sys.argv)):
    if sys.argv[i] == "-b":
        W1 = int(sys.argv[i+1])

if not W1:
    W1 = domain2.protect_sw_func.max_threads_per_block
W2 = 1
W3 = 1

get_kernel_function_info(domain2.protect_sw_func, 
        W1,W2, W3)


domain2.protect_sw_func(
    numpy.int32( domain2.number_of_elements),
    numpy.float64( domain2.minimum_allowed_height),
    numpy.float64( domain2.maximum_allowed_speed),
    numpy.float64( domain2.epsilon),
    drv.In( domain2.quantities['stage'].centroid_values),
    drv.In( domain2.quantities['elevation'].centroid_values),
    drv.InOut( domain2.quantities['xmomentum'].centroid_values),
    drv.InOut( domain2.quantities['ymomentum'].centroid_values),
    block = (W1, W2, W3),
    grid = ( (N + W1*W2*W3 -1)/(W1*W2*W3), 1)
)

cnt_x = 0
cnt_y = 0

x1 = domain1.quantities['xmomentum'].centroid_values
y1 = domain1.quantities['ymomentum'].centroid_values
x2 = domain2.quantities['xmomentum'].centroid_values
y2 = domain2.quantities['ymomentum'].centroid_values

for i in range(N):
    if x1[i] != x2[i]:
        cnt_x +=1
        if cnt_x < 5:
            print i, x1[i], x2[i]
    if y1[i] != y2[i]:
        cnt_y += 1
        if cnt_y < 5:
            print i, y1[i], y2[i]

print cnt_x, cnt_y
