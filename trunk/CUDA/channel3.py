"""Simple water flow example using ANUGA

Water flowing down a channel with more complex topography
"""

#---------------------------------------------------------------------------
# Import necessary modules
#---------------------------------------------------------------------------
import sys
import anuga
from anuga_cuda import GPU_domain

finaltime = 16.0

def generate_channel3_domain(gpu=True):
    #-----------------------------------------------------------------------
    # Setup computational domain
    #-----------------------------------------------------------------------
    length = 40.
    width = 5.
    dx = dy = .1           # Resolution: Length of subdivisions on both axes
    
    points, vertices, boundary = anuga.rectangular_cross(int(length/dx),
                                        int(width/dy), len1=length, len2=width)
    if gpu:
        domain = GPU_domain(points, vertices, boundary )
        if '-gpu' in sys.argv:
            domain.using_gpu = True
            print " --> Enable GPU version"
        for i in range(len(sys.argv)):
            if sys.argv[i] == '-fs':
                global finaltime
                finaltime = float(sys.argv[i+1])
                print " --> Finaltime is reset as %f" % finaltime

    else:
        domain = anuga.Domain(points, vertices, boundary)
    domain.set_name('channel3')                  # Output name
    #print domain.statistics()
    
    #-----------------------------------------------------------------------
    # Setup initial conditions
    #-----------------------------------------------------------------------
    def topography(x,y):
        """Complex topography defined by a function of vectors x and y."""
    
        z = -x/10
    
        N = len(x)
        for i in range(N):
            # Step
            if 10 < x[i] < 12:
                z[i] += 0.4 - 0.05*y[i]
    
            # Constriction
            if 27 < x[i] < 29 and y[i] > 3:
                z[i] += 2
    
            # Pole
            if (x[i] - 34)**2 + (y[i] - 2)**2 < 0.4**2:
                z[i] += 2
    
        return z
    
    domain.set_quantity('elevation', topography)    # elevation is a function
    domain.set_quantity('friction', 0.01)                  # Constant friction
    domain.set_quantity('stage', expression='elevation') # Dry initial condition
    
    #-----------------------------------------------------------------------
    # Setup boundary conditions
    #-----------------------------------------------------------------------
    Bi = anuga.Dirichlet_boundary([0.4, 0, 0])          # Inflow
    Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
    Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # Outflow
    
    domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})
    return domain

def evolve_channel3_domain(domain):
    #-----------------------------------------------------------------------
    # Evolve system through time
    #-----------------------------------------------------------------------
    global finaltime
    for t in domain.evolve(yieldstep=0.1, finaltime=finaltime):
        print domain.timestepping_statistics()
    
        #if domain.get_quantity('stage').\
        #       get_values(interpolation_points=[[10, 2.5]]) > 0:
        #    print 'Stage > 0: Changing to outflow boundary'
        #    domain.set_boundary({'right': Bo})
            
if __name__ == '__main__':
    domain = generate_channel3_domain()
    evolve_channel3_domain(domain)
