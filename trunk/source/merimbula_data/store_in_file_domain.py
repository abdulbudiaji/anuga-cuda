
def store_domain_in_file(file_name = "Merimbula_domain", domain=None):
    if(domain == None ):
        from generate_domain import domain_create
        domain = domain_create()

    
