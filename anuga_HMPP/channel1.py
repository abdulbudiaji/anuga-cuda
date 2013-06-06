"""Simple water flow example using ANUGA

Water flowing down a channel
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

# Import standard shallow water domain and standard boundaries.
import sys
import anuga
from hmpp_domain import HMPP_domain

finaltime= 1
yieldstep= 0.2

def generate_channel1_domain(gpu=True):
    #--------------------------------------------------------------------------
    # Setup computational domain
    #--------------------------------------------------------------------------

    #points, vertices, boundary = anuga.rectangular_cross(10, 1,
    #                                               len1=10.0, len2=5.0) # Mesh
    points, vertices, boundary = anuga.rectangular_cross(90, 90,
                                                   len1=0.1, len2=0.1) # Mesh
    
    if gpu:
        domain = HMPP_domain(points, vertices, boundary)  # Create domain
    else:
        domain = anuga.Domain(points, vertices, boundary)  # Create domain
    domain.set_name('channel1')                  # Output name
    
    #--------------------------------------------------------------------------
    # Setup initial conditions
    #--------------------------------------------------------------------------
    def topography(x, y):
        return -x/10                             # linear bed slope
    
    domain.set_quantity('elevation', topography) # Use function for elevation
    domain.set_quantity('friction', 0.01)        # Constant friction 
    domain.set_quantity('stage',                 # Dry bed
                        expression='elevation')  
    
    #--------------------------------------------------------------------------
    # Setup boundary conditions
    #--------------------------------------------------------------------------
    Bi = anuga.Dirichlet_boundary([0.4, 0, 0])         # Inflow
    Br = anuga.Reflective_boundary(domain)             # Solid reflective wall
    
    domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})


    return domain




def evolve_channel1_domain( domain ):
    #-----------------------------------------------------------------------
    # Evolve system through time
    #-----------------------------------------------------------------------
    global finaltime
    for t in domain.evolve(yieldstep=0.2, finaltime=finaltime):
        print domain.flux_timestep
        print domain.timestepping_statistics()
    

if __name__ == '__main__':
    import numpy
    from anuga.config import epsilon
    from utility import *
    domain = generate_channel1_domain()
    print domain.number_of_elements
    #evolve_channel1_domain(domain)
    
    testing_domain = generate_channel1_domain(False)

    cf_method, flow_algorithm, ts_method = number_domain_method(domain)

    
    from hmpp_python_glue import *

    hmpp_distribute_to_vertices_and_edges(
            domain,
            numpy.float64( yieldstep ),
            numpy.float64( finaltime ),
            numpy.float64( 0.0 ),
            numpy.float64( epsilon ),
            False,
            
            numpy.int32( cf_method ),
            numpy.int32( flow_algorithm ),
            numpy.int32( ts_method ),
            0
            )

    #testing_domain.distribute_to_vertices_and_edges()
    testing_domain.extrapolate_second_order_sw()


    check_all(domain, testing_domain)
