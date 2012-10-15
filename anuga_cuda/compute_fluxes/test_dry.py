def is_dry(domain, k, i = 0):
    ki = k * 3 + i

    n = domain.neighbours[ki]

    ql = list()
    qr = list()

    ql.append( domain.quantities['stage'].edge_values[ki] )
    ql.append( domain.quantities['xmomentum'].edge_values[ki] )
    ql.append( domain.quantities['ymomentum'].edge_values[ki] )
    zl = domain.quantities['elevation'].edge_values[ki]

    if n < 0:
        m = - n -1
        qr.append( domain.quantities['stage'].boundary_values[m] )
        qr.append( domain.quantities['xmomentum'].boundary_values[m] )
        qr.append( domain.quantities['ymomentum'].boundary_values[m] )
        zr = zl
    else:
        m = domain.neighbour_edges[ki];
        nm = n * 3 + m

        qr.append( domain.quantities['stage'].edge_values[nm] )
        qr.append( domain.quantities['xmomentum'].edge_values[nm] )
        qr.append( domain.quantities['ymomentum'].edge_values[nm] )
        zr = domain.quantities['elevation'].edge_values[nm] )
        
    if fabs(ql[0]-zl) < domain.epsilon and fabs(qr[0] -zr) < domain.epsilon:
        if n >= 0:
            print "This is the dry cell"
        else:
            print "n<0"
