import numpy
from pycuda import driver as drv
from anuga_cuda.merimbula_data.channel1 import generate_domain

domain1 = generate_domain( gpu=False )
domain2 = generate_domain( gpu=True )

domain1.protect_against_infinitesimal_and_negative_heights()
domain2.protect_against_infinitesimal_and_negative_heights()

N = domain1.number_of_elements
W1 = 32
W2 = 1
W3 = 1

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
    for i in range(N):
        if (Q1.edge_values[i] != Q2.edge_values[i]).all():
            cnt_e += 1
            if cnt_e  < 5:
                print i, Q1.edge_values[i], Q2.edge_values[i]
        if (Q1.vertex_values[i] != Q2.vertex_values[i]).all():
            cnt_v += 1
            if cnt_v  < 5:
                print i, Q1.vertex_values[i], Q2.vertex_values[i]

    print "___________ %d, %d\n" % ( cnt_e, cnt_v)
