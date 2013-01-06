import numpy as num
from gpu_domain  import GPU_domain

"""test flux calculation (actual C implementation)

This one tests the constant case where only the pressure term
contributes to each edge and cancels out once the total flux has
been summed up.
"""

a = [0.0, 0.0]
b = [0.0, 2.0]
c = [2.0, 0.0]
d = [0.0, 4.0]
e = [2.0, 2.0]
f = [4.0, 0.0]

points = [a, b, c, d, e, f]
#              bac,     bce,     ecf,     dbe
vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

domain = GPU_domain(points, vertices)
domain.check_integrity()

# The constant case
domain.set_quantity('elevation', -1)
domain.set_quantity('stage', 1)

domain.evolve()
print "+++ fluxes: done +++"
domain.evolve()
print "+++ fluxes: done +++"
# Central triangle
#assert num.allclose(domain.get_quantity('stage').explicit_update[1], 0)

# The more general case
#def surface(x, y):
#    return -x/2

#domain.set_quantity('elevation', -10)
#domain.set_quantity('stage', surface)
#domain.set_quantity('xmomentum', 1)

#domain.compute_fluxes()

#print domain.get_quantity('stage').explicit_update
# FIXME (Ole): TODO the general case
#assert allclose(domain.get_quantity('stage').explicit_update[1], ...??)
