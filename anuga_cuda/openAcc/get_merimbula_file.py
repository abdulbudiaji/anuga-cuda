#!/usr/bin/env python

from anuga_cuda import generate_merimbula_domain

domain = generate_merimbula_domain()

fileHandle = open('merimbula_cf.dat', 'w')

fileHandle.write("%d\n" % domain.number_of_elements)
# # of boundary
fileHandle.write("%d\n"%domain.quantities['stage'].boundary_values.shape[0])
fileHandle.write("%d\n" % domain.optimise_dry_cells)

fileHandle.write("%lf\n" % domain.evolve_max_timestep)
fileHandle.write("%lf\n" % domain.g)
fileHandle.write("%lf\n" % domain.epsilon)
fileHandle.write("%lf\n" % (domain.H0 * domain.H0))
fileHandle.write("%lf\n" % domain.H0 * 10)


s = domain.quantities['stage']
b = domain.quantities['elevation']
x = domain.quantities['xmomentum']
y = domain.quantities['ymomentum']


sv = s.vertex_values    # N * 3
sb = s.boundary_values  # N2
se = s.edge_values      # N * 3
sc = s.centroid_values  # N
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
normals = domain.normals    # N * 6
edgelengths = domain.edgelengths    # N * 3
radii = domain.radii        # N
areas = domain.areas        # N
tri = domain.tri_full_flag  # N
vertex_coordinates = domain.vertex_coordinates  # N * 3 * 2


for i in range(domain.number_of_elements):
    fileHandle.write(
        "%lf %lf %lf\n" % (sv[i][0], sv[i][1], sv[i][2]))
    fileHandle.write(
        "%lf %lf %lf\n" % (se[i][0], se[i][1], se[i][2]))
    fileHandle.write( "%lf\n" % sc[i])
    fileHandle.write( "%lf\n" % su[i])

    fileHandle.write(
        "%lf %lf %lf\n" % (be[i][0], be[i][1], be[i][2]))
    fileHandle.write( "%lf\n" % bc[i])
    

    
    fileHandle.write(
        "%lf %lf %lf\n" % (xe[i][0], xe[i][1], xe[i][2]))
    fileHandle.write( "%lf\n" % xu[i])

    fileHandle.write(
        "%lf %lf %lf\n" % (ye[i][0], ye[i][1], ye[i][2]))
    fileHandle.write( "%lf\n" % yu[i])

    
    fileHandle.write(
        "%lf %lf %lf\n" % (n[i][0], n[i][1], n[i][2]))

    fileHandle.write(
        "%lf %lf %lf\n" % (ne[i][0], ne[i][1], ne[i][2]))
    
    fileHandle.write( "%lf %lf %lf %lf %lf %lf\n" %\
        (normals[i][0], normals[i][1], normals[i][2], 
        normals[i][3], normals[i][4], normals[i][5] ))

    fileHandle.write( "%lf %lf %lf\n" %\
        (edgelengths[i][0], edgelengths[i][1], edgelengths[i][2]))

    fileHandle.write( "%lf\n" % radii[i])
    fileHandle.write( "%lf\n" % areas[i])
    fileHandle.write( "%ld\n" % tri[i])


    fileHandle.write( "%lf %lf\n" % \
        (vertex_coordinates[i*3][0], vertex_coordinates[i*3][1]))
    fileHandle.write( "%lf %lf\n" % \
        (vertex_coordinates[i*3 + 1][0], vertex_coordinates[i*3 + 1][1]))
    fileHandle.write( "%lf %lf\n" % \
        (vertex_coordinates[i*3 + 2][0], vertex_coordinates[i*3 + 2][1]))

for i in range(sb.shape[0]):
    fileHandle.write( "%lf\n" % sb[i])
    fileHandle.write( "%lf\n" % xb[i])
    fileHandle.write( "%lf\n" % yb[i])
    

fileHandle.close()
