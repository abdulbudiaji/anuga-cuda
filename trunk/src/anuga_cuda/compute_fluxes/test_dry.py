from anuga_cuda.merimbula_data.generate_domain import *
domain = domain_create()

def is_dry(domain, k, i = 0):

    n = domain.neighbours[k][i]

    ql = list()
    qr = list()

    ql.append( domain.quantities['stage'].edge_values[k][i] )
    ql.append( domain.quantities['xmomentum'].edge_values[k][i] )
    ql.append( domain.quantities['ymomentum'].edge_values[k][i] )
    zl = domain.quantities['elevation'].edge_values[k][i]

    if n < 0:
        m = - n -1
        qr.append( domain.quantities['stage'].boundary_values[m] )
        qr.append( domain.quantities['xmomentum'].boundary_values[m] )
        qr.append( domain.quantities['ymomentum'].boundary_values[m] )
        zr = zl
    else:
        m = domain.neighbour_edges[k][i];
        

        qr.append( domain.quantities['stage'].edge_values[n][m] )
        qr.append( domain.quantities['xmomentum'].edge_values[n][m] )
        qr.append( domain.quantities['ymomentum'].edge_values[n][m] )
        zr = domain.quantities['elevation'].edge_values[n][m] 

    if abs(ql[0]-zl) < domain.epsilon and abs(qr[0] -zr) < domain.epsilon:
        if n >= 0:
            print "This is the dry cell"
        else:
            print "n<0"
