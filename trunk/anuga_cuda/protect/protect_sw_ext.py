import numpy
from pycuda import driver as drv
from anuga_cuda import generate_merimbula_domain

domain1 = generate_merimbula_domain()
domain2 = generate_merimbula_domain(gpu=True)

for i in domain1.evolve(yieldstep = 50, finaltime = 50):
    pass
domain1.protect_against_infinitesimal_and_negative_heights()

for i in domain2.evolve(yieldstep = 50, finaltime = 50):
    pass
N = domain2.number_of_elements
W1 = 32
W2 = 1
W3 = 1
domain2.equip_kernel_functions()

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

print x1, y1
print x2, y2

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
