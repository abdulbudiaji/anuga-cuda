#!/usr/bin/env python


def get_domain_for_compute_fluxes( domain, fileHandle):
    # Number of elements
    fileHandle.write("%ld " % domain.number_of_elements)
    # Number of boundaries
    fileHandle.write("%ld "%domain.quantities['stage'].boundary_values.shape[0])
    # Number of boundary tags
    fileHandle.write("%ld " % domain.boundary_map.__len__())
    
    
    fileHandle.write("%lf " % domain.epsilon)
    fileHandle.write("%lf " % domain.H0)
    fileHandle.write("%lf " % domain.g)
    fileHandle.write("%d " % domain.optimise_dry_cells)
    fileHandle.write("%lf " % domain.evolve_max_timestep)
    fileHandle.write("%lf " % domain.extrapolate_velocity_secnod_order)
    fileHandle.write("%lf " % domain.minimum_allowed_height)
    
    
    fileHandle.write("%lf " % domain.beta_w )
    fileHandle.write("%lf " % domain.beta_w_dry )
    fileHandle.write("%lf " % domain.beta_uh )
    fileHandle.write("%lf " % domain.beta_uh_dry )
    fileHandle.write("%lf " % domain.beta_vh )
    fileHandle.write("%lf " % domain.beta_vh_dry )
    
    fileHandle.write("%s\n" % domain.compute_fluxes_method)
    
    s = domain.quantities['stage']
    b = domain.quantities['elevation']
    xm = domain.quantities['xmomentum']
    ym = domain.quantities['ymomentum']
    h = domain.quantities['height']
    f = domain.quantities['friction']
    xv = domain.quantities['xvelocity']
    yv = domain.quantities['yvelocity']
    
    
    
    se = s.edge_values      # N * 3
    sc = s.centroid_values  # N
    sv = s.vertex_values    # N * 3
    sb = s.boundary_values  # N2
    su = s.explicit_update  # N
    
    be = b.edge_values      # N * 3
    bc = b.centroid_values  # N
    
    xb = x.boundary_values  # N2
    xe = x.edge_values      # N * 3
    xu = x.explicit_update  # N
    
    yb = y.boundary_values  # N2
    ye = y.edge_values      # N * 3
    yu = y.explicit_update  # N
    
    
    n = domain.neighbours       # N * 3
    ne = domain.neighbour_edges # N * 3
    sn = domain.surrogate_neighbours
    normals = domain.normals    # N * 6
    edgelengths = domain.edgelengths    # N * 3
    radii = domain.radii        # N
    areas = domain.areas        # N
    tri = domain.tri_full_flag  # N
    maxspeed = domain.max_speed
    
    
    vertex_coordinates = domain.vertex_coordinates  # N * 3 * 2
    edge_coordinates = domain.vertex_coordinates  # N * 3 * 2
    centroid_coordinates = domain.vertex_coordinates  # N * 3 * 2
    
    nb = domain.number_of_boundaries
    
    
    for i in range(domain.number_of_elements):
        fileHandle.write(
            "%ld %ld %ld\n" % (n[i][0], n[i][1], n[i][2]))
    
        fileHandle.write(
            "%ld %ld %ld\n" % (ne[i][0], ne[i][1], ne[i][2]))
        
        fileHandle.write(
            "%ld %ld %ld\n" % (sn[i][0], sn[i][1], sn[i][2]))
        
        fileHandle.write( "%lf %lf %lf %lf %lf %lf\n" %\
            (normals[i][0], normals[i][1], normals[i][2], 
            normals[i][3], normals[i][4], normals[i][5] ))
    
        fileHandle.write( "%lf %lf %lf\n" %\
            (edgelengths[i][0], edgelengths[i][1], edgelengths[i][2]))
    
        fileHandle.write( "%lf\n" % radii[i])
        fileHandle.write( "%lf\n" % areas[i])
        fileHandle.write( "%ld\n" % tri[i])
        #fileHandle.write( "%lf\n" % maxspeed[i])
    
    
    
        fileHandle.write( "%lf %lf %lf %lf %lf \n" % \
            (vertex_coordinates[i*3][0], vertex_coordinates[i*3][1],
            vertex_coordinates[i*3 + 1][0], vertex_coordinates[i*3 + 1][1],
            vertex_coordinates[i*3 + 2][0], vertex_coordinates[i*3 + 2][1]))
    
        fileHandle.write( "%lf %lf %lf %lf %lf \n" % \
            (edge_coordinates[i*3][0], edge_coordinates[i*3][1],
            edge_coordinates[i*3 + 1][0], edge_coordinates[i*3 + 1][1],
            edge_coordinates[i*3 + 2][0], edge_coordinates[i*3 + 2][1]))
    
        fileHandle.write( "%lf %lf %lf %lf %lf \n" % \
            (centroid_coordinates[i*3][0], centroid_coordinates[i*3][1],
            centroid_coordinates[i*3 + 1][0], centroid_coordinates[i*3 + 1][1],
            centroid_coordinates[i*3 + 2][0], centroid_coordinates[i*3 + 2][1]))
    
    
        fileHandle.write("%ld \n" % nb[i])
        
        # Edge values
        fileHandle.write(
            "%lf %lf %lf\n" % (se[i][0], se[i][1], se[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (xme[i][0], xme[i][1], xme[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (yme[i][0], yme[i][1], yme[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (be[i][0], be[i][1], be[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (he[i][0], he[i][1], he[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (xve[i][0], xve[i][1], xve[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (yve[i][0], yve[i][1], yve[i][2]))
        
    
        # Centroid values
        fileHandle.write(
            "%lf %lf %lf\n" % (sc[i][0], sc[i][1], sc[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (xmc[i][0], xmc[i][1], xmc[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (ymc[i][0], ymc[i][1], ymc[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (bc[i][0], bc[i][1], bc[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (fc[i][0], fc[i][1], fc[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (hc[i][0], hc[i][1], hc[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (xvc[i][0], xvc[i][1], xvc[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (yvc[i][0], yvc[i][1], yvc[i][2]))
        
    
        # Vertex values
        fileHandle.write(
            "%lf %lf %lf\n" % (sv[i][0], sv[i][1], sv[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (xmv[i][0], xmv[i][1], xmv[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (ymv[i][0], ymv[i][1], ymv[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (bv[i][0], bv[i][1], bv[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (hv[i][0], hv[i][1], hv[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (xvv[i][0], xvv[i][1], xvv[i][2]))
        fileHandle.write(
            "%lf %lf %lf\n" % (yvv[i][0], yvv[i][1], yvv[i][2]))
        
    
        # Explicit update
        fileHandle.write( "%lf\n" % su[i])
        fileHandle.write( "%lf\n" % xmu[i])
        fileHandle.write( "%lf\n" % ymu[i])
    
        # Semi-implicit update
        fileHandle.write( "%lf\n" % sse[i])
        fileHandle.write( "%lf\n" % xmse[i])
        fileHandle.write( "%lf\n" % ymse[i])
    
        
        # Gradient values
        fileHandle.write( "%lf %lf %lf %lf %lf %lf \n" % (
            sxg[i], xmxg[i], ymxg[i], hxg[i], xvxg[i], yvxg[i]))
        fileHandle.write( "%lf %lf %lf %lf %lf %lf \n" % (
            syg[i], xmyg[i], ymyg[i], hyg[i], xvyg[i], yvyg[i]))
    
        
        # Some others
        fileHandle.write( "%lf %lf %d \n" % ( minbed[i], maxbed[i], cntwet[i]))
    
    
    
    for i in range(sb.shape[0]):
        fileHandle.write( "%lf %lf %lf %lf %lf %lf %lf \n" % (
            sb[i], xmb[i], ymb[i], bb[i], hb[i], xvb[i], yvb[i]))
        fileHandle.write( "%ld %ld \n" % (bcell[i], bedges[i]))
    
    
    for tag in domain.tag_boundary_cells:
        if self.boundary_map[tag] is None:
            fileHandle.write("0\n")
        
        ids = domain.tag_boundary_cells[tag]
        N = ids.__len__()
        bcell = domain.boundary_cells
        bedge = domain.boundary_edges
        for i in range(N):
            fileHandle.write("%ld %ld %ld\n" % (
                ids[i], bcell[ids[i]], bedge[ids[i]]))
        
def get_domain_for_ext2_lmt_byV( domain, fileHandle):
    # Number of elements
    fileHandle.write("%ld\n %lf\n" % (domain.number_of_elements, 
            domain.quantities['stage'].beta))
    
    n = domain.neighbours       # N * 3
    sn = domain.surrogate_neighbours # N * 3
    nb = domain.number_of_boundaries # N
    
    centroid_coordinates = domain.centroid_coordinates  # N * 2
    vertex_coordinates = domain.vertex_coordinates  # N * 3 * 2
    
    
    
    s = domain.quantities['stage']
    
    
    
    se = s.edge_values      # N * 3
    sc = s.centroid_values  # N
    sv = s.vertex_values    # N * 3
    
    
    
    for i in range(domain.number_of_elements):
        fileHandle.write( "%lf %lf\n" % \
            (centroid_coordinates[i][0], centroid_coordinates[i][1]))
        
        fileHandle.write( "%lf %lf %lf %lf %lf %lf\n" % \
            (vertex_coordinates[i*3][0], vertex_coordinates[i*3][1],
            vertex_coordinates[i*3 + 1][0], vertex_coordinates[i*3 + 1][1],
            vertex_coordinates[i*3 + 2][0], vertex_coordinates[i*3 + 2][1]))

        fileHandle.write("%d\n" % nb[i])

        fileHandle.write(
            "%ld %ld %ld\n" % (sn[i][0], sn[i][1], sn[i][2]))
        
        fileHandle.write(
            "%ld %ld %ld\n" % (n[i][0], n[i][1], n[i][2]))
    
    
        
        # Centroid values
        fileHandle.write(
            "%lf\n" % sc[i])
        
    
        # Vertex values
        fileHandle.write(
            "%lf %lf %lf\n" % (sv[i][0], sv[i][1], sv[i][2]))


        # Edge values
        fileHandle.write(
            "%lf %lf %lf\n" % (se[i][0], se[i][1], se[i][2]))
        
    
       
    
    
    
    
        

if __name__ == "__main__":
    from merimbula import generate_merimbula_domain
    
    domain = generate_merimbula_domain()
    
    fileHandle = open('merimbula_all.dat', 'w')
    
    get_domain_for_ext2_lmt_byV(domain, fileHandle)
    
    fileHandle.close()
