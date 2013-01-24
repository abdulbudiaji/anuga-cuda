
import numpy
from pycuda import driver as drv
from anuga_cuda.merimbula_data.channel1 import generate_domain
from anuga_cuda.merimbula_data.generate_domain import domain_create

#domain1 = generate_domain( gpu=False )
#domain2 = generate_domain( gpu=True )
domain1 = domain_create( gpu=False )
domain2 = domain_create( gpu=True )

domain1.distribute_to_vertices_and_edges()
domain2.distribute_to_vertices_and_edges()


    
N = domain1.quantities['elevation'].boundary.shape[0]
W1 = 32
W2 = 1
W3 = 1

domain1.update_centroids_of_velocities_and_height_func()

domain2.quantities['elevation'].set_boundary_values_from_edges()

domain2.interpolate_from_vertices_to_edges_func(
    numpy.int32( N ),
    drv.In( domain2.quantities['stage'].centroid_values ),
    drv.In( domain2.quantities['xmomentum'].centroid_values ),
    drv.In( domain2.quantities['ymomentum'].centroid_values ),
    drv.InOut( domain2.quantities['height'].centroid_values ),
    drv.InOut( domain2.quantities['elevation'].centroid_values ),
    drv.InOut( domain2.quantities['xvelocity'].centroid_values ),
    drv.In( domain2.quantities['yvelocity'].centroid_values ),

    drv.In( domain2.quantities['stage'].boundary_values ),
    drv.In( domain2.quantities['xmomentum'].boundary_values ),
    drv.In( domain2.quantities['ymomentum'].boundary_values ),
    drv.InOut( domain2.quantities['height'].boundary_values ),
    drv.In( domain2.quantities['elevation'].boundary_values ),
    drv.In( domain2.quantities['xvelocity'].boundary_values ),
    drv.In( domain2.quantities['yvelocity'].boundary_values ),

    block = (W1, W2, W3),
    grid = ( (N + W1*W2*W3 -1) / (W1*W2*W3), 1)
    )

print numpy.allclose(Q1.edge_values, Q2.edge_values)
