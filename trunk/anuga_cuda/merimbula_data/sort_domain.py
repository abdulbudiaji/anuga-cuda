#!/usr/bin/env python

def swap_domain(domain, i, j, k):
    neighbours = domain.neighbours[i]
    n_edges = domain.neighbour_edges[i]
    s_neighbours = domain.surrogate_neighbours[i]
    normals = domain.normals[i]
    edgelengths = domain.edgelengths[i]

    vertex_coordinates = domain.vertex_coordinates
    edge_coordinates = domain.edge_coordinates


    stage_edge_values = domain.quantities['stage'].edge_values[i]
    xmom_edge_values = domain.quantities['xmomentum'].edge_values[i]
    ymom_edge_values = domain.quantities['ymomentum'].edge_values[i]
    bed_edge_values = domain.quantities['elevation'].edge_values[i]

    stage_vertex_values = domain.quantities['stage'].vertex_values[i]
    xmom_vertex_values = domain.quantities['xmomentum'].vertex_values[i]
    ymom_vertex_values = domain.quantities['ymomentum'].vertex_values[i]
    bed_vertex_values = domain.quantities['elevation'].vertex_values[i]
    

    neighbours[j], neighbours[k] = neighbours[k], neighbours[j]

    n_edges[j], n_edges[k] = n_edges[k], n_edges[j]

    s_neighbours[j], s_neighbours[k] = s_neighbours[k], s_neighbours[j]

    normals[j*2], normals[k*2] = normals[k*2], normals[j*2]

    normals[j*2+1], normals[k*2+1] = normals[k*2+1], normals[j*2+1]

    edgelengths[j], edgelengths[k] = edgelengths[k], edgelengths[j]

    vertex_coordinates[i*3+j][0], vertex_coordinates[i*3+k][0] = vertex_coordinates[i*3+k][0], vertex_coordinates[i*3+j][0] 
    vertex_coordinates[i*3+j][1], vertex_coordinates[i*3+k][1] = vertex_coordinates[i*3+k][1], vertex_coordinates[i*3+j][1] 

    edge_coordinates[i*3+j][0], edge_coordinates[i*3+k][0] = edge_coordinates[i*3+k][0], edge_coordinates[i*3+j][0] 
    edge_coordinates[i*3+j][1], edge_coordinates[i*3+k][1] = edge_coordinates[i*3+k][1], edge_coordinates[i*3+j][1] 

    stage_edge_values[j], stage_edge_values[k] = stage_edge_values[k], stage_edge_values[j]
    xmom_edge_values[j], xmom_edge_values[k] = xmom_edge_values[k], xmom_edge_values[j]
    ymom_edge_values[j], ymom_edge_values[k] = ymom_edge_values[k], ymom_edge_values[j]
    bed_edge_values[j], bed_edge_values[k] = bed_edge_values[k], bed_edge_values[j]

    stage_vertex_values[j], stage_vertex_values[k] = stage_vertex_values[k], stage_vertex_values[j]
    xmom_vertex_values[j], xmom_vertex_values[k] = xmom_vertex_values[k], xmom_vertex_values[j]
    ymom_vertex_values[j], ymom_vertex_values[k] = ymom_vertex_values[k], ymom_vertex_values[j]
    bed_vertex_values[j], bed_vertex_values[k] = bed_vertex_values[k], bed_vertex_values[j]



def sort_domain(domain):
    for i in range(domain.number_of_elements):
        
        neighbours = domain.neighbours[i]

        if neighbours[2]>=0 and neighbours[2] < i and \
                (neighbours[1]<0 or neighbours[2]<neighbours[1]):
            swap_domain(domain, i, 2, 1)

        if neighbours[1]>=0 and neighbours[1]<i and \
                (neighbours[0]<0 or neighbours[1]<neighbours[0]):
            swap_domain(domain, i, 0, 1)

            # changes
        if neighbours[2]>=0 and neighbours[2]<i and \
               (neighbours[1]<0 or neighbours[2]<neighbours[1]):
            swap_domain(domain, i, 2, 1)
    


if __name__ == "__main__":
    from anuga_cuda.merimbula_data.generate_domain import domain_create

    from anuga_cuda.compute_fluxes.compute_fluxes import spe_bubble_sort

    domain1 = domain_create()
    domain2 = domain_create()
    
    sort_domain(domain2)

    for k in range(domain1.number_of_elements):
        b = [0,1,2]
        spe_bubble_sort(b, domain1.neighbours[k], k)
        for i in range(3):
            if domain1.neighbours[k][b[i]] != domain2.neighbours[k][i]:
                print "###### tri: %ld, edge: %d - %d " % (k, b[i], i)
                continue
            
            if domain1.neighbour_edges[k][ b[i] ] != domain2.neighbour_edges[k][i]:
                print "tri: %ld, edge: %d - %d @ neighbour_edges " % (k, b[i], i)

            if domain1.surrogate_neighbours[k][ b[i] ] != domain2.surrogate_neighbours[k][i]:
                print "tri: %ld, edge: %d - %d @ surrogate_neighbours " % (k, b[i], i)
                
            if domain1.normals[k][ b[i]*2] != domain2.normals[k][i*2] or\
                     domain1.normals[k][ b[i]*2+1] != domain2.normals[k][i*2+1]:
                print "tri: %ld, edge: %ld - %d @ normals " % (k, b[i], i)

            if domain1.edgelengths[k][ b[i] ] != domain2.edgelengths[k][i]:
                print "tri: %ld, edge: %ld - %d @ edgelengths " % (k, b[i], i)


            if domain1.vertex_coordinates[k*3+b[i]][0] != domain2.vertex_coordinates[k*3+i][0] or\
                     domain1.vertex_coordinates[k*3+b[i]][1] != domain2.vertex_coordinates[k*3+i][1]:
                print "tri: %ld, edge: %ld - %d @ vertex_coordinates " % (k, b[i], i)

            if domain1.edge_coordinates[k*3+b[i]][0] != domain2.edge_coordinates[k*3+i][0] or\
                     domain1.edge_coordinates[k*3+b[i]][1] != domain2.edge_coordinates[k*3+i][1]:
                print "tri: %ld, edge: %ld - %d @ edge_coordinates " % (k, b[i], i)


            if domain1.quantities['stage'].edge_values[k][ b[i] ] != \
                    domain2.quantities['stage'].edge_values[k][i]:
                print "tri: %ld, edge: %ld - %d @ stage_edge_values " % (k, b[i], i)
    
            if domain1.quantities['xmomentum'].edge_values[k][ b[i] ] != \
                    domain2.quantities['xmomentum'].edge_values[k][i]:
                print "tri: %ld, edge: %ld - %d @ xmom_edge_values " % (k, b[i], i)

            if domain1.quantities['ymomentum'].edge_values[k][ b[i] ] != \
                    domain2.quantities['ymomentum'].edge_values[k][i]:
                print "tri: %ld, edge: %ld - %d @ ymom_edge_values " % (k, b[i], i)

            if domain1.quantities['elevation'].edge_values[k][ b[i] ] != \
                    domain2.quantities['elevation'].edge_values[k][i]:
                print "tri: %ld, edge: %ld - %d @ bed_edge_values " % (k, b[i], i)


            if domain1.quantities['stage'].vertex_values[k][ b[i] ] != \
                    domain2.quantities['stage'].vertex_values[k][i]:
                print "tri: %ld, vertex: %ld - %d @ stage_vertex_values " % (k, b[i], i)
    
            if domain1.quantities['xmomentum'].vertex_values[k][ b[i] ] != \
                    domain2.quantities['xmomentum'].vertex_values[k][i]:
                print "tri: %ld, vertex: %ld - %d @ xmom_vertex_values " % (k, b[i], i)

            if domain1.quantities['ymomentum'].vertex_values[k][ b[i] ] != \
                    domain2.quantities['ymomentum'].vertex_values[k][i]:
                print "tri: %ld, vertex: %ld - %d @ ymom_vertex_values " % (k, b[i], i)

            if domain1.quantities['elevation'].vertex_values[k][ b[i] ] != \
                    domain2.quantities['elevation'].vertex_values[k][i]:
                print "tri: %ld, vertex: %ld - %d @ bed_vertex_values " % (k, b[i], i)


