#!/usr/bin/env python

def swap_domain(domain, i, j, k):
    neighbours = domain.neighbours[i]
    
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

        if neighbours[2]>=0 and neighbours[2]<i and \
               (neighbours[1]<0 or neighbours[2]<neighbours[1]):
            swap_domain(domain, i, 2, 1)
    
    
    for i in range(domain.number_of_elements):
        n1, n2, n3 = domain.neighbours[i]
        
        for index in range(3):
            if domain.neighbours[n1][index] == i:
                domain.neighbour_edges[n1][index] = 0
            if domain.neighbours[n2][index] == i:
                domain.neighbour_edges[n2][index] = 1
            if domain.neighbours[n3][index] == i:
                domain.neighbour_edges[n3][index] = 2
                
                    


def rearrange_domain(domain):
    import copy 

    new_dom = copy.deepcopy(domain)

    sort_domain(new_dom)

    temp_dom = copy.deepcopy(new_dom)

    N = new_dom.number_of_elements

    neighbours = new_dom.neighbours[i]
    neighbour_edges = new_dom.neighbour_edges[i]
    surrogate_neighbours = new_dom.surrogate_neighbours[i]
    normals = new_dom.normals[i]
    edgelengths = new_dom.edgelengths[i]

    vertex_coordinates = new_dom.vertex_coordinates
    edge_coordinates = new_dom.edge_coordinates


    stage_edge_values = new_dom.quantities['stage'].edge_values[i]
    xmom_edge_values = new_dom.quantities['xmomentum'].edge_values[i]
    ymom_edge_values = new_dom.quantities['ymomentum'].edge_values[i]
    bed_edge_values = new_dom.quantities['elevation'].edge_values[i]

    stage_vertex_values = new_dom.quantities['stage'].vertex_values[i]
    xmom_vertex_values = new_dom.quantities['xmomentum'].vertex_values[i]
    ymom_vertex_values = new_dom.quantities['ymomentum'].vertex_values[i]
    bed_vertex_values = new_dom.quantities['elevation'].vertex_values[i]


    for i in range(N):
        neighbours[i/3][i%3] = temp_dom.neighbours[i][0]
        neighbours[(i+N)/3][(i+N)%3] = temp_dom.neighbours[i][1]
        neighbours[(i+2*N)/3][(i+2*N)%3] = temp_dom.neighbours[i][2]

        neighbour_edges[i/3][i%3] = temp_dom.neighbour_edges[i][0]
        neighbour_edges[(i+N)/3][(i+N)%3] = temp_dom.neighbour_edges[i][1]
        neighbour_edges[(i+2*N)/3][(i+2*N)%3] = temp_dom.neighbour_edges[i][2]

        surrogate_neighbours[i/3][i%3] = temp_dom.surrogate_neighbours[i][0]
        surrogate_neighbours[(i+N)/3][(i+N)%3] = temp_dom.surrogate_neighbours[i][1]
        surrogate_neighbours[(i+2*N)/3][(i+2*N)%3] = temp_dom.surrogate_neighbours[i][2]

        normals[i/3][i%3] = temp_dom.normals[i][0]
        normals[(i+N)/3][(i+N)%3] = temp_dom.normals[i][1]
        normals[(i+2*N)/3][(i+2*N)%3] = temp_dom.normals[i][2]
        normals[(i+3*N)/3][(i+3*N)%3] = temp_dom.normals[i][3]
        normals[(i+4*N)/3][(i+4*N)%3] = temp_dom.normals[i][4]
        normals[(i+5*N)/3][(i+5*N)%3] = temp_dom.normals[i][5]

        edgelengths[i/3][i%3] = temp_dom.edgelengths[i][0]
        edgelengths[(i+N)/3][(i+N)%3] = temp_dom.edgelengths[i][1]
        edgelengths[(i+2*N)/3][(i+2*N)%3] = temp_dom.edgelengths[i][2]

        
        vertex_coordinates[i/2][i%2] = temp_dom.vertex_coordinates[i][0]
        vertex_coordinates[(i+N)/2][(i+N)%2] = temp_dom.vertex_coordinates[i][1]
        vertex_coordinates[(i+N*2)/2][(i+N*2)%2] = temp_dom.vertex_coordinates[i][2]
        vertex_coordinates[(i+N*3)/2][(i+N*3)%2] = temp_dom.vertex_coordinates[i][3]
        vertex_coordinates[(i+N*4)/2][(i+N*4)%2] = temp_dom.vertex_coordinates[i][4]
        vertex_coordinates[(i+N*5)/2][(i+N*5)%2] = temp_dom.vertex_coordinates[i][5]

        
        edge_coordinates[i/2][i%2] = temp_dom.edge_coordinates[i][0]
        edge_coordinates[(i+N)/2][(i+N)%2] = temp_dom.edge_coordinates[i][1]
        edge_coordinates[(i+N*2)/2][(i+N*2)%2] = temp_dom.edge_coordinates[i][2]
        edge_coordinates[(i+N*3)/2][(i+N*3)%2] = temp_dom.edge_coordinates[i][3]
        edge_coordinates[(i+N*4)/2][(i+N*4)%2] = temp_dom.edge_coordinates[i][4]
        edge_coordinates[(i+N*5)/2][(i+N*5)%2] = temp_dom.edge_coordinates[i][5]


        centroid_coordinates[i/2][i%2] = temp_dom.centroid_coordinates[i][0]
        centroid_coordinates[(i+N)/2][(i+N)%2] = temp_dom.centroid_coordinates[i][1]

        stage_edge_values[i/3][i%3] = temp_dom.stage_edge_values[i][0]
        stage_edge_values[(i+N)/3][(i+N)%3] = temp_dom.stage_edge_values[i][1]
        stage_edge_values[(i+2*N)/3][(i+2*N)%3] = temp_dom.stage_edge_values[i][2]

        xmom_edge_values[i/3][i%3] = temp_dom.xmom_edge_values[i][0]
        xmom_edge_values[(i+N)/3][(i+N)%3] = temp_dom.xmom_edge_values[i][1]
        xmom_edge_values[(i+2*N)/3][(i+2*N)%3] = temp_dom.xmom_edge_values[i][2]

        ymom_edge_values[i/3][i%3] = temp_dom.ymom_edge_values[i][0]
        ymom_edge_values[(i+N)/3][(i+N)%3] = temp_dom.ymom_edge_values[i][1]
        ymom_edge_values[(i+2*N)/3][(i+2*N)%3] = temp_dom.ymom_edge_values[i][2]

        bed_edge_values[i/3][i%3] = temp_dom.bed_edge_values[i][0]
        bed_edge_values[(i+N)/3][(i+N)%3] = temp_dom.bed_edge_values[i][1]
        bed_edge_values[(i+2*N)/3][(i+2*N)%3] = temp_dom.bed_edge_values[i][2]



        stage_vertex_values[i/3][i%3] = temp_dom.stage_vertex_values[i][0]
        stage_vertex_values[(i+N)/3][(i+N)%3] = temp_dom.stage_vertex_values[i][1]
        stage_vertex_values[(i+2*N)/3][(i+2*N)%3] = temp_dom.stage_vertex_values[i][2]

        xmom_vertex_values[i/3][i%3] = temp_dom.xmom_vertex_values[i][0]
        xmom_vertex_values[(i+N)/3][(i+N)%3] = temp_dom.xmom_vertex_values[i][1]
        xmom_vertex_values[(i+2*N)/3][(i+2*N)%3] = temp_dom.xmom_vertex_values[i][2]

        ymom_vertex_values[i/3][i%3] = temp_dom.ymom_vertex_values[i][0]
        ymom_vertex_values[(i+N)/3][(i+N)%3] = temp_dom.ymom_vertex_values[i][1]
        ymom_vertex_values[(i+2*N)/3][(i+2*N)%3] = temp_dom.ymom_vertex_values[i][2]

        bed_vertex_values[i/3][i%3] = temp_dom.bed_vertex_values[i][0]
        bed_vertex_values[(i+N)/3][(i+N)%3] = temp_dom.bed_vertex_values[i][1]
        bed_vertex_values[(i+2*N)/3][(i+2*N)%3] = temp_dom.bed_vertex_values[i][2]


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


