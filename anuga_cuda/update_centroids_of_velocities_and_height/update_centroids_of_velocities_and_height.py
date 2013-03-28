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
    domain1 = generate_merimbula_domain(False)
    domain2 = generate_merimbula_domain(True)

if using_rearranged_domain:
    domain2 = rearrange_domain(domain2)
    sort_domain(domain1)

domain2.equip_kernel_functions()

N = domain2.number_of_elements
import sys
W1 = 0
for i in range( len(sys.argv)):
    if sys.argv[i] == "-b":
        W1 = int(sys.argv[i+1])

if not W1:
    W1 = domain2.update_centroids_of_velocities_and_height_func.max_threads_per_block
W2 = 1
W3 = 1

get_kernel_function_info(domain2.update_centroids_of_velocities_and_height_func, 
        W1,W2, W3)


print "Testing input"
print "  stage centroid_values all close? ", numpy.allclose( 
        domain1.quantities['stage'].centroid_values,
        domain2.quantities['stage'].centroid_values)
print "  xmomentum centroid_values all close? ", numpy.allclose( 
        domain1.quantities['xmomentum'].centroid_values,
        domain2.quantities['xmomentum'].centroid_values)
print "  ymomentum centroid_values all close? ", numpy.allclose( 
        domain1.quantities['ymomentum'].centroid_values,
        domain2.quantities['ymomentum'].centroid_values)
print "  height centroid_values all close? ", numpy.allclose( 
        domain1.quantities['height'].centroid_values,
        domain2.quantities['height'].centroid_values)
print "  elevation centroid_values all close? ", numpy.allclose( 
        domain1.quantities['elevation'].centroid_values,
        domain2.quantities['elevation'].centroid_values)
print "  xvelocity centroid_values all close? ", numpy.allclose( 
        domain1.quantities['xvelocity'].centroid_values,
        domain2.quantities['xvelocity'].centroid_values)
print "  yvelocity centroid_values all close? ", numpy.allclose( 
        domain1.quantities['yvelocity'].centroid_values,
        domain2.quantities['yvelocity'].centroid_values)

print "  stage boundary_values all close? ", numpy.allclose( 
        domain1.quantities['stage'].boundary_values,
        domain2.quantities['stage'].boundary_values)
print "  xmomentum boundary_values all close? ", numpy.allclose( 
        domain1.quantities['xmomentum'].boundary_values,
        domain2.quantities['xmomentum'].boundary_values)
print "  ymomentum boundary_values all close? ", numpy.allclose( 
        domain1.quantities['ymomentum'].boundary_values,
        domain2.quantities['ymomentum'].boundary_values)
print "  height boundary_values all close? ", numpy.allclose( 
        domain1.quantities['height'].boundary_values,
        domain2.quantities['height'].boundary_values)
print "  elevation boundary_values all close? ", numpy.allclose( 
        domain1.quantities['elevation'].boundary_values,
        domain2.quantities['elevation'].boundary_values)
print "  xvelocity boundary_values all close? ", numpy.allclose( 
        domain1.quantities['xvelocity'].boundary_values,
        domain2.quantities['xvelocity'].boundary_values)
print "  yvelocity boundary_values all close? ", numpy.allclose( 
        domain1.quantities['yvelocity'].boundary_values,
        domain2.quantities['yvelocity'].boundary_values)



domain1.update_centroids_of_velocities_and_height()

domain2.quantities['elevation'].set_boundary_values_from_edges()
domain2.update_centroids_of_velocities_and_height_func(
        numpy.int32(N),
        numpy.int32(domain2.quantities['stage'].boundary_values.shape[0]),

        drv.In( domain2.quantities['stage'].centroid_values ),
        drv.In( domain2.quantities['xmomentum'].centroid_values ),
        drv.In( domain2.quantities['ymomentum'].centroid_values ),
        drv.InOut( domain2.quantities['height'].centroid_values ),
        drv.In( domain2.quantities['elevation'].centroid_values ),
        drv.InOut( domain2.quantities['xvelocity'].centroid_values ),
        drv.InOut( domain2.quantities['yvelocity'].centroid_values ),

        drv.In( domain2.quantities['stage'].boundary_values ),
        drv.In( domain2.quantities['xmomentum'].boundary_values ),
        drv.In( domain2.quantities['ymomentum'].boundary_values ),
        drv.InOut( domain2.quantities['height'].boundary_values ),
        drv.In( domain2.quantities['elevation'].boundary_values ),
        drv.InOut( domain2.quantities['xvelocity'].boundary_values ),
        drv.InOut( domain2.quantities['yvelocity'].boundary_values ),
        block = (W1, W2, W3),
        grid=((N+W1*W2*W3-1)/(W1*W2*W3),1)
        )

hc1 = domain1.quantities['height'].centroid_values
hb1 = domain1.quantities['height'].boundary_values
uc1 = domain1.quantities['xmomentum'].centroid_values
ub1 = domain1.quantities['xmomentum'].boundary_values
vc1 = domain1.quantities['ymomentum'].centroid_values
vb1 = domain1.quantities['ymomentum'].boundary_values

hc2 = domain2.quantities['height'].centroid_values
hb2 = domain2.quantities['height'].boundary_values
uc2 = domain2.quantities['xmomentum'].centroid_values
ub2 = domain2.quantities['xmomentum'].boundary_values
vc2 = domain2.quantities['ymomentum'].centroid_values
vb2 = domain2.quantities['ymomentum'].boundary_values
wb2 = domain2.quantities['stage'].boundary_values
zb2 = domain2.quantities['elevation'].boundary_values


res = []
res.append( numpy.allclose( hc1, hc2) )
res.append( numpy.allclose( hb1, hb2) )
res.append( numpy.allclose( uc1, uc2) )
res.append( numpy.allclose( ub1, ub2) )
res.append( numpy.allclose( vc1, vc2) )
res.append( numpy.allclose( vb1, vb2) )

print "\n\nTesting output"
print "  height centroid_values all closed? " ,res[0]
print "  height boundary_values all closed? " ,res[1]
print "  xmom centroid_values all closed? " ,res[2]
print "  xmom boundary_values all closed? " ,res[3]
print "  ymom centroid_values all closed? " ,res[4]
print "  ymom boundary_values all closed? " ,res[5]

if not res.count(True) == 6:
    cnt_hb = 0
    res = []
    for i in range( domain1.quantities['stage'].boundary_values.shape[0] ):
        if hb1[i] != hb2[i]:
            cnt_hb += 1
            res.append(hb2[i])
            if cnt_hb <10:
                print "# %d, %f, %f, wb=%f, zb=%f" % \
                    (i, hb1[i], hb2[i], wb2[i], zb2[i])
    print " Number of diff: is %d " % cnt_hb
    #print res
