
import numpy
from pycuda import driver as drv
from anuga_cuda import generate_merimbula_domain
from anuga_cuda import generate_channel3_domain
from anuga_cuda import generate_cairns_domain, get_kernel_function_info

using_tsunami_domain = False

if using_tsunami_domain:
    domain1 = generate_cairns_domain(False)
    domain2 = generate_cairns_domain(True)
else:
    domain1 = generate_merimbula_domain(False)
    domain2 = generate_merimbula_domain(True)

domain2.equip_kernel_functions()


domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
domain1.evolve(yieldstep = 50, finaltime = 500)
    

domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)
domain2.evolve(yieldstep = 50, finaltime = 500)


domain1.protect_against_infinitesimal_and_negative_heights()
domain2.protect_against_infinitesimal_and_negative_heights()


if domain1.optimised_gradient_limiter:
    if domain1._order_ == 1:
        for name in domain1.conserved_quantities:
            Q1 = domain1.quantities[name]
            Q1.extrapolate_first_order()
            
            Q2 = domain2.quantities[name]
            Q2.extrapolate_first_order()
    
    
    elif domain1._order_ == 2:
        domain1.extrapolate_second_order_sw()
        domain2.extrapolate_second_order_sw()
else:
    for name in domain1.conserved_quantities:
        Q1 = domain1.quantities[name]
        Q2 = domain1.quantities[name]
        if domain1._order_ == 1:
            Q1.extrapolate_first_order()
            Q2.extrapolate_first_order()
        elif domain1._order_ == 2:
            Q1.extrapolate_second_order_and_limit_by_vertex()
            Q2.extrapolate_second_order_and_limit_by_vertex()
        
N = domain1.number_of_elements
import sys
W1 = 0
for i in range( len(sys.argv)):
    if sys.argv[i] == "-b":
        W1 = int(sys.argv[i+1])

if not W1:
    W1 = domain2.balance_deep_and_shallow_func.max_threads_per_block
print W1
W2 = 1
W3 = 1

get_kernel_function_info(domain2.balance_deep_and_shallow_func, W1,W2, W3)


domain1.balance_deep_and_shallow()
domain2.balance_deep_and_shallow_func(
    numpy.int32(N),
    numpy.float64( domain2.H0 ),
    numpy.float64( domain2.alpha_balance ),
    numpy.int32( domain2.tight_slope_limiters ),
    numpy.int32( domain2.use_centroid_velocities),

    drv.In( domain2.quantities['stage'].centroid_values ),
    drv.In( domain2.quantities['elevation'].centroid_values ),
    drv.InOut( domain2.quantities['stage'].vertex_values ),
    drv.In( domain2.quantities['elevation'].vertex_values ),


    drv.In( domain2.quantities['xmomentum'].centroid_values ),
    drv.In( domain2.quantities['ymomentum'].centroid_values ),
    drv.InOut( domain2.quantities['xmomentum'].vertex_values ),
    drv.InOut( domain2.quantities['ymomentum'].vertex_values ),

    block = (W1, W2, W3),
    grid = ( (N+ W1*W2*W3 - 1)/(W1*W2*W3), 1)
    )


cnt_sv = 0
cnt_ev = 0
cnt_xv = 0
cnt_yv = 0

sv1 = domain1.quantities['stage'].vertex_values
sv2 = domain2.quantities['stage'].vertex_values

ev1 = domain1.quantities['elevation'].vertex_values
ev2 = domain2.quantities['elevation'].vertex_values

xv1 = domain1.quantities['xmomentum'].vertex_values
xv2 = domain2.quantities['xmomentum'].vertex_values

yv1 = domain1.quantities['ymomentum'].vertex_values
yv2 = domain2.quantities['ymomentum'].vertex_values

print numpy.allclose( sv1, sv2 )
print numpy.allclose( ev1, ev2 )
print numpy.allclose( xv1, xv2 )
print numpy.allclose( yv1, yv2 )

cnt_sc = 0
cnt_ec = 0
cnt_xc = 0
cnt_yc = 0

sc1 = domain1.quantities['stage'].centroid_values
sc2 = domain2.quantities['stage'].centroid_values

ec1 = domain1.quantities['elevation'].centroid_values
ec2 = domain2.quantities['elevation'].centroid_values

xc1 = domain1.quantities['xmomentum'].centroid_values
xc2 = domain2.quantities['xmomentum'].centroid_values

yc1 = domain1.quantities['ymomentum'].centroid_values
yc2 = domain2.quantities['ymomentum'].centroid_values

print numpy.allclose( sc1, sc2)
print numpy.allclose( ec1, ec2)
print numpy.allclose( xc1, xc2)
print numpy.allclose( yc1, yc2)

for i in range(N):
    if (sv1[i] != sv2[i]).all():
        cnt_sv += 1
        if cnt_sv <= 5:
            print i, sv1[i], sv2[i]
    if (ev1[i] != ev2[i]).all():
        cnt_ev += 1
    if (xv1[i] != xv2[i]).all():
        cnt_xv += 1
    if (yv1[i] != yv2[i]).all():
        cnt_yv += 1

print cnt_sv, cnt_ev, cnt_xv, cnt_yv


for i in range(N):
    if (sc1[i] != sc2[i]).all():
        cnt_sc += 1
    if (ec1[i] != ec2[i]).all():
        cnt_ec += 1
    if (xc1[i] != xc2[i]).all():
        cnt_xc += 1
    if (yc1[i] != yc2[i]).all():
        cnt_yc += 1

print cnt_sc, cnt_ec, cnt_xc, cnt_yc

