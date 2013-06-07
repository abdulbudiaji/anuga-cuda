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
                    


def sort_domain_check(domain2):
    from anuga_cuda.merimbula_data.generate_domain import domain_create

    domain1 = domain_create()

    for k in range(domain1.number_of_elements):
        b = [0,1,2]
        spe_bubble_sort(b, domain1.neighbours[k], k)
        for i in range(3):
            if domain1.neighbours[k][b[i]] != domain2.neighbours[k][i]:
                print "###### tri: %ld, edge: %d - %d " % (k, b[i], i)
                continue
            
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



def rearrange_domain(domain, sort_flag=True, spare_domain=None):
    import copy 

    if spare_domain == None:
        new_dom = copy.deepcopy(domain)
	
        if sort_flag:
            sort_domain(new_dom)
            temp_dom = copy.deepcopy(new_dom)
        else:
            temp_dom = domain
    else:
        new_dom = domain
        temp_dom = spare_domain

    N = new_dom.number_of_elements

    neighbours = new_dom.neighbours
    neighbour_edges = new_dom.neighbour_edges
    surrogate_neighbours = new_dom.surrogate_neighbours
    normals = new_dom.normals
    edgelengths = new_dom.edgelengths

    vertex_coordinates = new_dom.vertex_coordinates
    edge_coordinates = new_dom.edge_coordinates
    centroid_coordinates = new_dom.centroid_coordinates


    stage_edge_values = new_dom.quantities['stage'].edge_values
    xmom_edge_values = new_dom.quantities['xmomentum'].edge_values
    ymom_edge_values = new_dom.quantities['ymomentum'].edge_values
    bed_edge_values = new_dom.quantities['elevation'].edge_values

    stage_vertex_values = new_dom.quantities['stage'].vertex_values
    xmom_vertex_values = new_dom.quantities['xmomentum'].vertex_values
    ymom_vertex_values = new_dom.quantities['ymomentum'].vertex_values
    bed_vertex_values = new_dom.quantities['elevation'].vertex_values


    for i in range(N):
        neighbours[i/3][i%3] = temp_dom.neighbours[i][0]
        neighbours[(i+N)/3][(i+N)%3] = temp_dom.neighbours[i][1]
        neighbours[(i+2*N)/3][(i+2*N)%3] = temp_dom.neighbours[i][2]

        neighbour_edges[i/3][i%3] = temp_dom.neighbour_edges[i][0]
        neighbour_edges[(i+N)/3][(i+N)%3] = temp_dom.neighbour_edges[i][1]
        neighbour_edges[(i+2*N)/3][(i+2*N)%3] = temp_dom.neighbour_edges[i][2]

        surrogate_neighbours[i/3][i%3] = \
										temp_dom.surrogate_neighbours[i][0]
        surrogate_neighbours[(i+N)/3][(i+N)%3] = \
										temp_dom.surrogate_neighbours[i][1]
        surrogate_neighbours[(i+2*N)/3][(i+2*N)%3] = \
										temp_dom.surrogate_neighbours[i][2]


        normals[i/6][i%6] = temp_dom.normals[i][0]
        normals[(i+N)/6][(i+N)%6] = temp_dom.normals[i][1]
        normals[(i+2*N)/6][(i+2*N)%6] = temp_dom.normals[i][2]
        normals[(i+3*N)/6][(i+3*N)%6] = temp_dom.normals[i][3]
        normals[(i+4*N)/6][(i+4*N)%6] = temp_dom.normals[i][4]
        normals[(i+5*N)/6][(i+5*N)%6] = temp_dom.normals[i][5]


        edgelengths[i/3][i%3] = temp_dom.edgelengths[i][0]
        edgelengths[(i+N)/3][(i+N)%3] = temp_dom.edgelengths[i][1]
        edgelengths[(i+2*N)/3][(i+2*N)%3] = temp_dom.edgelengths[i][2]

        
        #vertex_coordinates[i/2][i%2] = \
		#							temp_dom.vertex_coordinates[i*3][0]
        #vertex_coordinates[(i+N)/2][(i+N)%2] = \
		#							temp_dom.vertex_coordinates[i*3][1]
        #vertex_coordinates[(i+N*2)/2][(i+N*2)%2] = \
		#							temp_dom.vertex_coordinates[i*3+1][0]
        #vertex_coordinates[(i+N*3)/2][(i+N*3)%2] = \
		#							temp_dom.vertex_coordinates[i*3+1][1]
        #vertex_coordinates[(i+N*4)/2][(i+N*4)%2] = \
		#							temp_dom.vertex_coordinates[i*3+2][0]
        #vertex_coordinates[(i+N*5)/2][(i+N*5)%2] = \
        #							temp_dom.vertex_coordinates[i*3+2][1]

        vertex_coordinates.flat[i] = \
									temp_dom.vertex_coordinates[i*3][0]
        vertex_coordinates.flat[i+N] = \
									temp_dom.vertex_coordinates[i*3][1]
        vertex_coordinates.flat[i+N*2] = \
									temp_dom.vertex_coordinates[i*3+1][0]
        vertex_coordinates.flat[i+N*3] = \
									temp_dom.vertex_coordinates[i*3+1][1]
        vertex_coordinates.flat[i+N*4] = \
									temp_dom.vertex_coordinates[i*3+2][0]
        vertex_coordinates.flat[i+N*5] = \
        							temp_dom.vertex_coordinates[i*3+2][1]

        
        edge_coordinates[i/2][i%2] = temp_dom.edge_coordinates[i*3][0]
        edge_coordinates[(i+N)/2][(i+N)%2] = temp_dom.edge_coordinates[i*3][1]
        edge_coordinates[(i+N*2)/2][(i+N*2)%2] = temp_dom.edge_coordinates[i*3+1][0]
        edge_coordinates[(i+N*3)/2][(i+N*3)%2] = temp_dom.edge_coordinates[i*3+1][1]
        edge_coordinates[(i+N*4)/2][(i+N*4)%2] = temp_dom.edge_coordinates[i*3+2][0]
        edge_coordinates[(i+N*5)/2][(i+N*5)%2] = temp_dom.edge_coordinates[i*3+2][1]


        centroid_coordinates[i/2][i%2] = temp_dom.centroid_coordinates[i][0]
        centroid_coordinates[(i+N)/2][(i+N)%2] = temp_dom.centroid_coordinates[i][1]

        stage_edge_values[i/3][i%3] = \
							temp_dom.quantities['stage'].edge_values[i][0]
        stage_edge_values[(i+N)/3][(i+N)%3] = \
							temp_dom.quantities['stage'].edge_values[i][1]
        stage_edge_values[(i+2*N)/3][(i+2*N)%3] = \
							temp_dom.quantities['stage'].edge_values[i][2]


        xmom_edge_values[i/3][i%3] = \
						temp_dom.quantities['xmomentum'].edge_values[i][0]
        xmom_edge_values[(i+N)/3][(i+N)%3] = \
						temp_dom.quantities['xmomentum'].edge_values[i][1]
        xmom_edge_values[(i+2*N)/3][(i+2*N)%3] = \
                        temp_dom.quantities['xmomentum'].edge_values[i][2]


        ymom_edge_values[i/3][i%3] = \
							temp_dom.quantities['ymomentum'].edge_values[i][0]
        ymom_edge_values[(i+N)/3][(i+N)%3] = \
							temp_dom.quantities['ymomentum'].edge_values[i][1]
        ymom_edge_values[(i+2*N)/3][(i+2*N)%3] = \
							temp_dom.quantities['ymomentum'].edge_values[i][2]


        bed_edge_values[i/3][i%3] = \
							temp_dom.quantities['elevation'].edge_values[i][0]
        bed_edge_values[(i+N)/3][(i+N)%3] = \
							temp_dom.quantities['elevation'].edge_values[i][1]
        bed_edge_values[(i+2*N)/3][(i+2*N)%3] = \
							temp_dom.quantities['elevation'].edge_values[i][2]



        stage_vertex_values[i/3][i%3] = \
							temp_dom.quantities['stage'].vertex_values[i][0]
        stage_vertex_values[(i+N)/3][(i+N)%3] = \
							temp_dom.quantities['stage'].vertex_values[i][1]
        stage_vertex_values[(i+2*N)/3][(i+2*N)%3] = \
							temp_dom.quantities['stage'].vertex_values[i][2]


        xmom_vertex_values[i/3][i%3] = \
							temp_dom.quantities['xmomentum'].vertex_values[i][0]
        xmom_vertex_values[(i+N)/3][(i+N)%3] = \
							temp_dom.quantities['xmomentum'].vertex_values[i][1]
        xmom_vertex_values[(i+2*N)/3][(i+2*N)%3] = \
							temp_dom.quantities['xmomentum'].vertex_values[i][2]


        ymom_vertex_values[i/3][i%3] = \
							temp_dom.quantities['ymomentum'].vertex_values[i][0]
        ymom_vertex_values[(i+N)/3][(i+N)%3] = \
							temp_dom.quantities['ymomentum'].vertex_values[i][1]
        ymom_vertex_values[(i+2*N)/3][(i+2*N)%3] = \
							temp_dom.quantities['ymomentum'].vertex_values[i][2]


        bed_vertex_values[i/3][i%3] = \
							temp_dom.quantities['elevation'].vertex_values[i][0]
        bed_vertex_values[(i+N)/3][(i+N)%3] = \
							temp_dom.quantities['elevation'].vertex_values[i][1]
        bed_vertex_values[(i+2*N)/3][(i+2*N)%3] = \
							temp_dom.quantities['elevation'].vertex_values[i][2]

    return new_dom



def check_rearranged_array(a, b, configure=3):
    import numpy 
    counter = 0
    N = a.shape[0]
    if configure == 2:
        """ N * 2 """
        for i in range(N):
            if a[i][0] != b[i/2][i%2] or \
		    	a[i][1] != b[(i+N)/2][(i+N)%2] :
                counter +=1

    if configure == 3:
        """ N * 3 """
        for i in range(N):
            if not ( numpy.allclose( a[i][0], b[i/3][i%3] ) and \
	        	numpy.allclose( a[i][1], b[(i+N)/3][(i+N)%3] ) and \
                numpy.allclose( a[i][2], b[(i+N*2)/3][(i+N*2)%3] ) ):
                if counter < 5:
                    print i, a[i], b[i/3][i%3], b[(i+N)/3][(i+N)%3], \
                            b[(i+N*2)/3][(i+N*2)%3]

                counter += 1

    elif configure == 6:
        """ N * 6 """
        for i in range(N):
            if a[i][0] != b[i/6][i%6] or \
                a[i][1] != b[(i+N)/6][(i+N)%6] or \
	    	    a[i][2] != b[(i+N*2)/6][(i+N*2)%6] or \
	    	    a[i][3] != b[(i+N*3)/6][(i+N*3)%6] or \
	    	    a[i][4] != b[(i+N*4)/6][(i+N*4)%6] or \
	    	    a[i][5] != b[(i+N*5)/6][(i+N*5)%6]:

	    		counter +=1

    elif configure == 32:
        """ N*3 * 2"""
        N = N/3
        for i in range(N):
            if a[i*3][0] != b[i/2][i%2] or \
		    	a[i*3][1] != b[(i+N)/2][(i+N)%2] or \
		    	a[i*3 + 1][0] != b[(i+N*2)/2][(i+N*2)%2] or \
		    	a[i*3 + 1][1] != b[(i+N*3)/2][(i+N*3)%2] or \
		    	a[i*3 + 2][0] != b[(i+N*4)/2][(i+N*4)%2] or \
		    	a[i*3 + 2][1] != b[(i+N*5)/2][(i+N*5)%2] :
                
                counter +=1
    
    return counter


    
def rearrange_domain_check(domain1, domain2):
    neighbours = domain1.neighbours
    neighbour_edges = domain1.neighbour_edges
    surrogate_neighbours = domain1.surrogate_neighbours
    normals = domain1.normals
    edgelengths = domain1.edgelengths

    vertex_coordinates = domain1.vertex_coordinates
    edge_coordinates = domain1.edge_coordinates
    centroid_coordinates = domain1.centroid_coordinates


    stage_edge_values = domain1.quantities['stage'].edge_values
    xmom_edge_values = domain1.quantities['xmomentum'].edge_values
    ymom_edge_values = domain1.quantities['ymomentum'].edge_values
    bed_edge_values = domain1.quantities['elevation'].edge_values

    stage_vertex_values = domain1.quantities['stage'].vertex_values
    xmom_vertex_values = domain1.quantities['xmomentum'].vertex_values
    ymom_vertex_values = domain1.quantities['ymomentum'].vertex_values
    bed_vertex_values = domain1.quantities['elevation'].vertex_values


    vc = domain2.vertex_coordinates
    xc = domain2.edge_coordinates
    cc = domain2.centroid_coordinates


    se = domain2.quantities['stage'].edge_values
    xe = domain2.quantities['xmomentum'].edge_values
    ye = domain2.quantities['ymomentum'].edge_values
    be = domain2.quantities['elevation'].edge_values

    sv = domain2.quantities['stage'].vertex_values
    xv = domain2.quantities['xmomentum'].vertex_values
    yv = domain2.quantities['ymomentum'].vertex_values
    bv = domain2.quantities['elevation'].vertex_values



    N = domain1.number_of_elements

    cnt_nb = check_rearranged_array(
            neighbours, domain2.neighbours, 3)
    cnt_nbedge = check_rearranged_array(
            neighbour_edges, domain2.neighbour_edges, 3)
    cnt_surrnb = check_rearranged_array(
            surrogate_neighbours, domain2.surrogate_neighbours, 3)
    cnt_norm = check_rearranged_array(
            normals, domain2.normals, 6)
    cnt_el = check_rearranged_array(
            edgelengths, domain2.edgelengths, 3)
    cnt_vc = check_rearranged_array(
            vertex_coordinates, domain2.vertex_coordinates, 32)
    cnt_ec = check_rearranged_array(
            edge_coordinates, domain2.edge_coordinates, 32)
    
    edgeC = edge_coordinates
    ec_2 = domain2.edge_coordinates.flat
    for i in range(N):
        if edgeC[i*3][0] != ec_2[i] or edgeC[i*3][1] != ec_2[i+N] or \
            edgeC[i*3+1][0]!=ec_2[i+N*2] or edgeC[i*3+1][1]!=ec_2[i+N*3] or\
            edgeC[i*3+2][0]!=ec_2[i+N*4] or edgeC[i*3+2][1]!=ec_2[i+N*5]:
            print i

    cnt_cc = check_rearranged_array(
            centroid_coordinates, domain2.centroid_coordinates, 2)
    cnt_se = check_rearranged_array( stage_edge_values,
            domain2.quantities['stage'].edge_values, 3)
    cnt_xe = check_rearranged_array( xmom_edge_values,
            domain2.quantities['xmomentum'].edge_values, 3)
    cnt_ye = check_rearranged_array( ymom_edge_values,
            domain2.quantities['ymomentum'].edge_values, 3)
    cnt_be = check_rearranged_array( bed_edge_values,
            domain2.quantities['elevation'].edge_values, 3)
    cnt_sv = check_rearranged_array( stage_vertex_values,
            domain2.quantities['stage'].vertex_values, 3)
    cnt_xv = check_rearranged_array( xmom_vertex_values,
            domain2.quantities['xmomentum'].vertex_values, 3)
    cnt_yv = check_rearranged_array( ymom_vertex_values,
            domain2.quantities['ymomentum'].vertex_values, 3)
    cnt_bv = check_rearranged_array( bed_vertex_values,
            domain2.quantities["elevation"].vertex_values, 3)

    #for i in range(N):
    #    if neighbours[i][0] != domain2.neighbours[i/3][i%3] or \
	#		neighbours[i][1] != domain2.neighbours[(i+N)/3][(i+N)%3] or \
	#		neighbours[i][2] != domain2.neighbours[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, neighbours[i], domain2.neighbours[i/3][i%3], \
    #                    domain2.neighbours[(i+N)/3][(i+N)%3], \
    #                    domain2.neighbours[(i+N*2)/3][(i+N*2)%3]
    #        cnt_nb+=1

    #    if neighbour_edges[i][0] != domain2.neighbour_edges[i/3][i%3] or \
	#		neighbour_edges[i][1] != \
    #            domain2.neighbour_edges[(i+N)/3][(i+N)%3]or\
	#		neighbour_edges[i][2] != \
    #            domain2.neighbour_edges[(i+N*2)/3][(i+N*2)%3]:
	#		cnt_nbedge+=1
    #    
    #    if surrogate_neighbours[i][0] != \
    #            domain2.surrogate_neighbours[i/3][i%3]or\
	#		surrogate_neighbours[i][1] !=\
    #            domain2.surrogate_neighbours[(i+N)/3][(i+N)%3]or\
	#		surrogate_neighbours[i][2] !=\
    #            domain2.surrogate_neighbours[(i+N*2)/3][(i+N*2)%3]:
	#		cnt_surrnb +=1

    #    if normals[i][0] != domain2.normals[i/6][i%6] or \
	#		    normals[i][1] != domain2.normals[(i+N)/6][(i+N)%6] or \
	#		    normals[i][2] != domain2.normals[(i+N*2)/6][(i+N*2)%6] or \
	#		    normals[i][3] != domain2.normals[(i+N*3)/6][(i+N*3)%6] or \
	#		    normals[i][4] != domain2.normals[(i+N*4)/6][(i+N*4)%6] or \
	#		    normals[i][5] != domain2.normals[(i+N*5)/6][(i+N*5)%6]:
	#		cnt_norm +=1

    #    if edgelengths[i][0] != domain2.edgelengths[i/3][i%3] or \
	#		edgelengths[i][1] != domain2.edgelengths[(i+N)/3][(i+N)%3] or \
	#		edgelengths[i][2] != domain2.edgelengths[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, edgelengths[i], domain2.edgelengths[i/3][i%3], \
    #                    domain2.edgelengths[(i+N)/3][(i+N)%3], \
    #                    domain2.edgelengths[(i+N*2)/3][(i+N*2)%3]
    #        cnt_el +=1


    #    if vertex_coordinates[i*3][0] != \
    #            domain2.vertex_coordinates[i/2][i%2] or \
	#		vertex_coordinates[i*3][1] != \
    #            domain2.vertex_coordinates[(i+N)/2][(i+N)%2] or \
	#		vertex_coordinates[i*3 + 1][0] != \
    #            domain2.vertex_coordinates[(i+N*2)/2][(i+N*2)%2] or \
	#		vertex_coordinates[i*3 + 1][1] != \
    #            domain2.vertex_coordinates[(i+N*3)/2][(i+N*3)%2] or \
	#		vertex_coordinates[i*3 + 2][0] != \
    #            domain2.vertex_coordinates[(i+N*4)/2][(i+N*4)%2] or \
	#		vertex_coordinates[i*3 + 2][1] != \
    #            domain2.vertex_coordinates[(i+N*5)/2][(i+N*5)%2] :
    #        cnt_vc +=1


    #    if edge_coordinates[i*3][0] != \
    #            domain2.edge_coordinates[i/2][i%2] or \
	#		edge_coordinates[i*3][1] != \
    #            domain2.edge_coordinates[(i+N)/2][(i+N)%2] or \
	#		edge_coordinates[i*3 + 1][0] != \
    #            domain2.edge_coordinates[(i+N*2)/2][(i+N*2)%2] or \
	#		edge_coordinates[i*3 + 1][1] != \
    #            domain2.edge_coordinates[(i+N*3)/2][(i+N*3)%2] or \
	#		edge_coordinates[i*3 + 2][0] != \
    #            domain2.edge_coordinates[(i+N*4)/2][(i+N*4)%2] or \
	#		edge_coordinates[i*3 + 2][1] != \
    #            domain2.edge_coordinates[(i+N*5)/2][(i+N*5)%2] :
    #        cnt_ec +=1


    #    if centroid_coordinates[i][0] !=\
    #            domain2.centroid_coordinates[i/2][i%2] or \
	#		centroid_coordinates[i][1] != \
    #            domain2.centroid_coordinates[(i+N)/2][(i+N)%2] :
    #        cnt_cc +=1


    #    if stage_edge_values[i][0] != \
    #            domain2.quantities["stage"].edge_values[i/3][i%3] or \
	#	    stage_edge_values[i][1] != \
    #            domain2.quantities["stage"].edge_values[(i+N)/3][(i+N)%3] \
    #        or stage_edge_values[i][2] != \
    #        domain2.quantities["stage"].edge_values[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, stage_edge_values[i], \
    #                domain2.quantities["stage"].edge_values[i/3][i%3], \
    #                domain2.quantities["stage"].edge_values[(i+N)/3][(i+N)%3],\
    #                domain2.quantities["stage"].edge_values[(i+N*2)/3][(i+N*2)%3]
    #        cnt_se+=1

    #    if xmom_edge_values[i][0] != \
    #            domain2.quantities["xmomentum"].edge_values[i/3][i%3] or \
	#		xmom_edge_values[i][1] != \
    #            domain2.quantities["xmomentum"].edge_values[(i+N)/3][(i+N)%3] or \
	#		xmom_edge_values[i][2] != \
    #            domain2.quantities["xmomentum"].edge_values[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, xmom_edge_values[i], domain2.quantities["xmomentum"].edge_values[i/3][i%3], domain2.quantities["xmomentum"].edge_values[(i+N)/3][(i+N)%3], domain2.quantities["xmomentum"].edge_values[(i+N*2)/3][(i+N*2)%3]
    #        cnt_xe+=1

    #    if ymom_edge_values[i][0] != \
    #        domain2.quantities["ymomentum"].edge_values[i/3][i%3] or \
	#		ymom_edge_values[i][1] != \
    #            domain2.quantities["ymomentum"].edge_values[(i+N)/3][(i+N)%3] or \
	#		ymom_edge_values[i][2] != \
    #            domain2.quantities["ymomentum"].edge_values[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, ymom_edge_values[i], domain2.quantities["ymomentum"].edge_values[i/3][i%3], domain2.quantities["ymomentum"].edge_values[(i+N)/3][(i+N)%3], domain2.quantities["ymomentum"].edge_values[(i+N*2)/3][(i+N*2)%3]
    #        cnt_ye+=1

    #    if bed_edge_values[i][0] != \
    #        domain2.quantities["elevation"].edge_values[i/3][i%3] or \
	#		bed_edge_values[i][1] != \
    #            domain2.quantities["elevation"].edge_values[(i+N)/3][(i+N)%3] or \
	#		bed_edge_values[i][2] != \
    #            domain2.quantities["elevation"].edge_values[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, bed_edge_values[i], domain2.quantities["elevation"].edge_values[i/3][i%3], domain2.quantities["elevation"].edge_values[(i+N)/3][(i+N)%3], domain2.quantities["elevation"].edge_values[(i+N*2)/3][(i+N*2)%3]
    #        cnt_be+=1

    #    # Vertex values

    #    if stage_vertex_values[i][0] != \
    #        domain2.quantities["stage"].vertex_values[i/3][i%3] or \
	#		stage_vertex_values[i][1] != \
    #            domain2.quantities["stage"].vertex_values[(i+N)/3][(i+N)%3] or \
	#		stage_vertex_values[i][2] != \
    #            domain2.quantities["stage"].vertex_values[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, stage_vertex_values[i], domain2.quantities["stage"].vertex_values[i/3][i%3], domain2.quantities["stage"].vertex_values[(i+N)/3][(i+N)%3], domain2.quantities["stage"].vertex_values[(i+N*2)/3][(i+N*2)%3]
    #        cnt_sv +=1

    #    if xmom_vertex_values[i][0] != \
    #        domain2.quantities["xmomentum"].vertex_values[i/3][i%3] or \
	#		xmom_vertex_values[i][1] != \
    #            domain2.quantities["xmomentum"].vertex_values[(i+N)/3][(i+N)%3] or \
	#		xmom_vertex_values[i][2] != \
    #            domain2.quantities["xmomentum"].vertex_values[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, xmom_vertex_values[i], domain2.quantities["xmomentum"].vertex_values[i/3][i%3], domain2.quantities["xmomentum"].vertex_values[(i+N)/3][(i+N)%3], domain2.quantities["xmomentum"].vertex_values[(i+N*2)/3][(i+N*2)%3]
    #        cnt_xv +=1

    #    if ymom_vertex_values[i][0] != \
    #        domain2.quantities["ymomentum"].vertex_values[i/3][i%3] or \
	#		ymom_vertex_values[i][1] != \
    #            domain2.quantities["ymomentum"].vertex_values[(i+N)/3][(i+N)%3] or \
	#		ymom_vertex_values[i][2] != \
    #            domain2.quantities["ymomentum"].vertex_values[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, ymom_vertex_values[i], domain2.quantities["ymomentum"].vertex_values[i/3][i%3], domain2.quantities["ymomentum"].vertex_values[(i+N)/3][(i+N)%3], domain2.quantities["ymomentum"].vertex_values[(i+N*2)/3][(i+N*2)%3]
    #        cnt_yv +=1

    #    if bed_vertex_values[i][0] != \
    #        domain2.quantities["elevation"].vertex_values[i/3][i%3] or \
	#		bed_vertex_values[i][1] != \
    #            domain2.quantities["elevation"].vertex_values[(i+N)/3][(i+N)%3] or \
	#		bed_vertex_values[i][2] != \
    #            domain2.quantities["elevation"].vertex_values[(i+N*2)/3][(i+N*2)%3]:
    #        if cnt_nb < 5:
    #            print i, bed_vertex_values[i], domain2.quantities["elevation"].vertex_values[i/3][i%3], domain2.quantities["elevation"].vertex_values[(i+N)/3][(i+N)%3], domain2.quantities["elevation"].vertex_values[(i+N*2)/3][(i+N*2)%3]
    #        cnt_bv +=1

    if cnt_nb or cnt_nbedge or cnt_surrnb or cnt_norm:
        print "nb=%d, nb_edge=%d surrnb=%d norm=%d" % \
            (cnt_nb, cnt_nbedge, cnt_surrnb, cnt_norm)
    if cnt_el or cnt_vc or cnt_ec or cnt_cc:
        print "elen=%d, vc=%d ec=%d cc=%d" % \
            (cnt_el, cnt_vc, cnt_ec, cnt_cc)
    if cnt_se or cnt_xe or cnt_ye or cnt_be:
        print "se=%d, xe=%d ye=%d be=%d" % \
            (cnt_se, cnt_xe, cnt_ye, cnt_be)
    if cnt_sv or cnt_xv or cnt_yv or cnt_bv:
        print "sv=%d, xv=%d yv=%d bv=%d" % \
            (cnt_sv, cnt_xv, cnt_yv, cnt_bv)




if __name__ == "__main__":
    from anuga_cuda import generate_merimbula_domain

    from anuga_cuda.compute_fluxes.compute_fluxes import spe_bubble_sort

    domain1 = generate_merimbula_domain()
    domain2 = generate_merimbula_domain()
    #sort_domain(domain2)
    #sort_domain_check(domain2)
	
    #sort_domain(domain1)
    print " **** rearrange_domain **** "
    domain2=rearrange_domain(domain2, False)
    print " **** check domain ****"
    rearrange_domain_check(domain1, domain2)
