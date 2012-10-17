def _rotate(q, n1, n2):
    q1 = q[1]
    q2 = q[2]

    q[1] = n1 * q1 + n2*q2
    q[2] = -n2 * q1 + n1*q2

def _compute_speed(uh, h, epsilon, h0, limiting_threshold):


    if h < limiting_threshold :
        
        if h < epsilon:
            h = 0.0
            u = 0.0
        else :
            u = uh / (h + h0 / h)
        
        uh = u * h
    else :
        u = uh / h
    return u



def _flux_function_central(q_left, q_right, z_left, z_right,
        n1, n2, epsilon, h0, limiting_threshold, g, edgeflux, max_speed):

    q_left_rotated = []
    q_right_rotated = []
    
    flux_right = []
    flux_left = []

    q_left_rotated.append( q_left[0] )
    q_right_rotated.append( q_right[0] )
    q_left_rotated.append( q_left[1] )
    q_right_rotated.append( q_right[1] )
    q_left_rotated.append( q_left[2] )
    q_right_rotated.append( q_right[2] )

    
    _rotate(q_left_rotated, n1, n2);
    _rotate(q_right_rotated, n1, n2);

    z = 0.5 * (z_left + z_right)

    w_left = q_left_rotated[0];
    h_left = w_left - z;
    uh_left = q_left_rotated[1];
    u_left = _compute_speed( uh_left, h_left,
            epsilon, h0, limiting_threshold);

    w_right = q_right_rotated[0];
    h_right = w_right - z;
    uh_right = q_right_rotated[1];
    u_right = _compute_speed( uh_right, h_right,
            epsilon, h0, limiting_threshold);

    
    vh_left = q_left_rotated[2];
    vh_right = q_right_rotated[2];

    
    
    
    _compute_speed(vh_left, h_left,
            epsilon, h0, limiting_threshold);
    _compute_speed(vh_right, h_right,
            epsilon, h0, limiting_threshold);

    
    soundspeed_left = sqrt(g * h_left);
    soundspeed_right = sqrt(g * h_right);

    s_max = max(u_left + soundspeed_left, u_right + soundspeed_right)

    if s_max < 0.0 :
        s_max = 0.0

    s_min = min(u_left - soundspeed_left, u_right - soundspeed_right)

    if s_min > 0.0:
        s_min = 0.0;

    flux_left[0] = u_left*h_left;
    flux_left[1] = u_left * uh_left + 0.5 * g * h_left*h_left;
    flux_left[2] = u_left*vh_left;

    flux_right[0] = u_right*h_right;
    flux_right[1] = u_right * uh_right + 0.5 * g * h_right*h_right;
    flux_right[2] = u_right*vh_right;

    
    denom = s_max - s_min;
    if denom < epsilon:
        memset(edgeflux, 0, 3 * sizeof (double));
        max_speed = 0.0;

    else :
        inverse_denominator = 1.0 / denom;
        for i in range(3) :
            edgeflux[i] = s_max * flux_left[i] - s_min * flux_right[i]
            edgeflux[i] += s_max * s_min * (q_right_rotated[i] - q_left_rotated[i])
            edgeflux[i] *= inverse_denominator
        

        
        max_speed = max(fabs(s_max), fabs(s_min));

        
        _rotate(edgeflux, n1, -n2);
    
    

def cf(domain, k , i):

    ql = []
    qr = []
    ql.append(domain.quantities['stage'].edge_values[k][i])
    ql.append(domain.quantities['xmomentum'].edge_values[k][i])
    ql.append(domain.quantities['ymomentum'].edge_values[k][i])
    zl = domain.quantities['elevation'].edge_values[k][i]

    n = domain.neighbours[k][i]

    if n < 0 :
        m = -n -1;
        qr.append(domain.quantities['stage'].boundary_values[m] )
        qr.append(domain.quantities['xmomentum'].boundary_values[m] )
        qr.append(domain.quantities['ymomentum'].boundary_values[m] )
        zr = zl
    else:
        m = domain.neighbour_edges[k][i]
        
        qr.append(domain.quantities['stage'].edge_values[n][m] )
        qr.append(domain.quantities['xmomentum'].edge_values[n][m] )
        qr.append(domain.quantities['ymomentum'].edge_values[n][m] )
        zr = domain.quantities['elevation'].edge_values[n][m]

    if( domain.optimise_dry_cells):
        if fabs(ql[0] - zl) < domain.epsilon and fabs(qr[0]-zr) < domain.epsilon:
            print "dry!!!!!!!"
            return

    _flux_function_central(ql, qr, zl, zr,
                    domain.normals[ki2], domain.normals[ki2 + 1],
                    domain.epsilon, domain.H0 * domain.H0, limiting_threshold, domain.g,
                    edgeflux, max_speed)



