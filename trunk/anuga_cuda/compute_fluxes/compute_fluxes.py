
compute_fluxes_central_call = 1
compute_fluxes_central_structure_call = 1

def compute_fluxes_central(domain):
    N = domain.quantities['stage'].edge_values.shap[0]
    W = 16
    compute_fluxes_central_call += 1
    
    elements = numpy.random.randn(6)
    elements = elements.astype(numpy.float64)
    elements[0] = domain.timestep
    elements[1] = domain.epsilon
    elements[2] = domain.H0
    elements[3] = domain.g
    elements[4] = domain.optimise_dry_cells
    elements[5] = compute_fluxes_central_call
    
    memset(domain.quantities['stage'].explicit_update, 0 , domain.quantities['stage'].explicit_update.shape[0] )
    memset(domain.quantities['xmomentum'].explicit_update, 0 , domain.quantities['xmomentum'].explicit_update.shape[0] )
    memset(domain.quantities['ymomentum'].explicit_update, 0 , domain.quantities['ymomentum'].explicit_update.shape[0] )
     
    mod = SourceModule("""
        __global__ void compute_fluxes_central(
                double *elements,
                long* neighbours,
                long* neighbour_edges,
                double* normals,
                double* edgelengths,
                double* radii,
                double* areas,
                long* tri_full_flag,
                double* stage_edge_values,
                double* xmom_edge_values,
                double* ymom_edge_values,
                double* bed_edge_values,
                double* stage_boundary_values,
                double* xmom_boundary_values,
                double* ymom_boundary_values,
                double* stage_explicit_update,
                double* xmom_explicit_update,
                double* ymom_explicit_update,
                long* already_computed_flux,
                double* max_speed_array)
        {
            int k = threadIdx.x + threadIdx.y;
            
            double max_speed, length, inv_area, zl, zr;
            double elements[2] = elements[2]*elements[2]; // This ensures a good balance when h approaches elements[2].

            double limiting_threshold = 10 * elements[2]; // Avoid applying limiter below this

            int i, m, n;
            int ki, nm = 0, ki2; // Index shorthands
            
            double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes
            
            for (i = 0; i < 3; i++) 
            {
                ki = k * 3 + i; // Linear index to edge i of triangle k
    
                if (already_computed_flux[ki] == call) {
                    // We've already computed the flux across this edge
                    continue;
                }
    
                // Get left hand side values from triangle k, edge i
                ql[0] = stage_edge_values[ki];
                ql[1] = xmom_edge_values[ki];
                ql[2] = ymom_edge_values[ki];
                zl = bed_edge_values[ki];
    
                // Get right hand side values either from neighbouring triangle
                // or from boundary array (Quantities at neighbour on nearest face).
                n = neighbours[ki];
                if (n < 0) {
                    // Neighbour is a boundary condition
                    m = -n - 1; // Convert negative flag to boundary index
    
                    qr[0] = stage_boundary_values[m];
                    qr[1] = xmom_boundary_values[m];
                    qr[2] = ymom_boundary_values[m];
                    zr = zl; // Extend bed elevation to boundary
                } else {
                    // Neighbour is a real triangle
                    m = neighbour_edges[ki];
                    nm = n * 3 + m; // Linear index (triangle n, edge m)
    
                    qr[0] = stage_edge_values[nm];
                    qr[1] = xmom_edge_values[nm];
                    qr[2] = ymom_edge_values[nm];
                    zr = bed_edge_values[nm];
                }
    
                // Now we have values for this edge - both from left and right side.
    
                if (elements[4]) {
                    // Check if flux calculation is necessary across this edge
                    // This check will exclude dry cells.
                    // This will also optimise cases where zl != zr as
                    // long as both are dry
    
                    if (fabs(ql[0] - zl) < elements[1] &&
                            fabs(qr[0] - zr) < elements[1]) {
                        // Cell boundary is dry
    
                        already_computed_flux[ki] = call; // #k Done
                        if (n >= 0) {
                            already_computed_flux[nm] = call; // #n Done
                        }
    
                        max_speed = 0.0;
                        continue;
                    }
                }
    
    
                if (fabs(zl - zr) > 1.0e-10) {
                    report_python_error(AT, "Discontinuous Elevation");
                    return 0.0;
                }
    
                // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
                ki2 = 2 * ki; //k*6 + i*2
    
                // Edge flux computation (triangle k, edge i)
                _flux_function_central(ql, qr, zl, zr,
                        normals[ki2], normals[ki2 + 1],
                        elements[1], elements[2], limiting_threshold, elements[3],
                        edgeflux, &max_speed);
    
    
                // Multiply edgeflux by edgelength
                length = edgelengths[ki];
                edgeflux[0] *= length;
                edgeflux[1] *= length;
                edgeflux[2] *= length;
    
    
                // Update triangle k with flux from edge i
                stage_explicit_update[k] -= edgeflux[0];
                xmom_explicit_update[k] -= edgeflux[1];
                ymom_explicit_update[k] -= edgeflux[2];
    
                already_computed_flux[ki] = call; // #k Done
    
    
                // Update neighbour n with same flux but reversed sign
                if (n >= 0) {
                    stage_explicit_update[n] += edgeflux[0];
                    xmom_explicit_update[n] += edgeflux[1];
                    ymom_explicit_update[n] += edgeflux[2];
    
                    already_computed_flux[nm] = call; // #n Done
                }
    
                // Update timestep based on edge i and possibly neighbour n
                if (tri_full_flag[k] == 1) {
                    if (max_speed > elements[1]) {
                        // Apply CFL condition for triangles joining this edge (triangle k and triangle n)
    
                        // CFL for triangle k
                        elements[0] = min(elements[0], radii[k] / max_speed);
    
                        if (n >= 0) {
                            // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                            elements[0] = min(elements[0], radii[n] / max_speed);
                        }
    
                        // Ted Rigby's suggested less conservative version
                        //if (n>=0) {
                        //  elements[0] = min(elements[0], (radii[k]+radii[n])/max_speed);
                        //} else {
                        //  elements[0] = min(elements[0], radii[k]/max_speed);
                        // }
                    }
                }
    
            } // End edge i (and neighbour n)
    
    
            // Normalise triangle k by area and store for when all conserved
            // quantities get updated
            inv_area = 1.0 / areas[k];
            stage_explicit_update[k] *= inv_area;
            xmom_explicit_update[k] *= inv_area;
            ymom_explicit_update[k] *= inv_area;
    
    
            // Keep track of maximal speeds
            max_speed_array[k] = max_speed;
            
        }
    """)

def compute_fluxes_central_structure(domain):
    pass

if __name__ == '__main__':
    from anuga_cuda.merimbula_data.generate_domain import domain_create
    domain = domain_create()
    
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import numpy 