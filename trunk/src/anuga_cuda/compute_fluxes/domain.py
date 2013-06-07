import numpy
from anuga_cuda.merimbula_data.generate_domain import domain_create
domain1 = domain_create()
domain2 = domain_create()


from shallow_water_ext import compute_fluxes_ext_central_structure
domain1.flux_timestep = compute_fluxes_ext_central_structure(domain1)

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
from anuga_cuda import kernel_path as kp
compute_fluxes_mod = SourceModule(
    open(kp["compute_fluxes_dir"]+"compute_fluxes.cu").read() )
compute_fluxes_central_function = compute_fluxes_mod.get_function(
    "compute_fluxes_central_structure_CUDA")

from anuga_cuda.compute_fluxes.compute_fluxes import \
    compute_fluxes_central_structure_single as cs


N = domain2.number_of_elements
W1 = 32
W2 = 1
W3 = 1

domain = domain2
timestep_array = numpy.zeros( 
        (domain.number_of_elements) , 
        dtype=numpy.float64)
for i in range(N):
    timestep_array[i] = domain.evolve_max_timestep
compute_fluxes_central_function( 
        numpy.uint(domain.number_of_elements),
        numpy.float64(domain.g),
        numpy.float64(domain.epsilon),
        numpy.float64(domain.H0 * domain.H0),
        numpy.float64(domain.H0*10),
        numpy.uint(domain.optimise_dry_cells),

        cuda.InOut( timestep_array ), 
        cuda.In( domain.neighbours ), 
        cuda.In( domain.neighbour_edges ),
        cuda.In( domain.normals ), 
        cuda.In( domain.edgelengths ), 
        cuda.In( domain.radii ), 
        cuda.In( domain.areas ), 
        cuda.In( domain.tri_full_flag ),
        cuda.In( domain.quantities['stage'].edge_values ), 
        cuda.In( domain.quantities['xmomentum'].edge_values ), 
        cuda.In( domain.quantities['ymomentum'].edge_values ), 
        cuda.In( domain.quantities['elevation'].edge_values ), 
        cuda.In( domain.quantities['stage'].boundary_values ), 
        cuda.In( domain.quantities['xmomentum'].boundary_values ), 
        cuda.In( domain.quantities['ymomentum'].boundary_values ), 
        cuda.InOut( domain.quantities['stage'].explicit_update ), 
        cuda.InOut( domain.quantities['xmomentum'].explicit_update ), 
        cuda.InOut( domain.quantities['ymomentum'].explicit_update ),  
        cuda.InOut( domain.max_speed), 
        block = ( W1, W2, W3),
        grid = ((N+W1*W2*W3-1)/(W1*W2*W3), 1) 
        )

b = numpy.argsort(timestep_array)
domain.flux_timestep = timestep_array[b[0]] 

