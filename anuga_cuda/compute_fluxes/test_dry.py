def is_dry(domain, k, i = 0):
    ki = k * 3 + i

    n = domain.neighbours[ki]

    ql = []
    qr = []

    ql[0] = domain.quantities['stage'].edge_values[ki]
    ql[1] = domain.quantities['xmomentum'].edge_values[ki]
    ql[2] = domain.quantities['ymomentum'].edge_values[ki]
    zl = domain.quantities['elevation'].edge_values[ki]

    if n < 0:
        m = - n -1
        qr[0] = domain.quantities['stage'].boundary_values[m]
        qr[1] = domain.quantities['xmomentum'].boundary_values[m]
        qr[2] = domain.quantities['ymomentum'].boundary_values[m];
        zr = zl
    else:
        m = domain.neighbour_edges[ki];
        nm = n * 3 + m

        qr[0] = domain.quantities['stage'].edge_values[nm];
        qr[1] = domain.quantities['xmomentum'].edge_values[nm];
        qr[2] = domain.quantities['ymomentum'].edge_values[nm];
        zr = domain.quantities['elevation'].edge_values[nm];
        
    if fabs(ql[0]-zl) < domain.epsilon and fabs(qr[0] -zr) < domain.epsilon:
        if n >= 0:
            print "This is the dry cell"
        else:
            print "n<0"
