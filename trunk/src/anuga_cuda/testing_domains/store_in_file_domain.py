
def file_write(fd, content):
    fd.write(str(content) + '\n')

def file_write_no_enter(fd, content):
    fd.write(str(content))

def 

def store_domain_in_file(file_name = "Merimbula_domain.d", domain=None):
    if(domain == None ):
        from generate_domain import domain_create
        domain = domain_create()
    
    fd = open(file_name, 'w')

    file_write( domain.number_of_elements )

    file_write( domain.epsilon)

    file_write( domain.H0)

    file_write( domain.g )

    file_write( domain.optimise_dry_cells)

    file_write( domain.evolve_max_timestep)

    file_write( domain.minimum_allowed_height)

    file_write( domain.extrapolate_velocity_second_order)

    file_write( domain.beta_w )
    file_write( domain.beta_w_dry )
    file_write( domain.beta_uh )
    file_write( domain.beta_uh_dry )
    file_write( domain.beta_vh )
    file_write( domain.beta_vh_dry )

    for i,
        file_write( domain.neighours)
