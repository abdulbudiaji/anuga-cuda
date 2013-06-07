""" Common filenames and locations for topographic data, meshes and outputs.
    This file defines the parameters of the scenario you wish to run.
"""

import anuga
from anuga_cuda import cairns_dir
#------------------------------------------------------------------------------
# Define scenario as either slide or fixed_wave. Choose one.
#------------------------------------------------------------------------------
scenario = 'fixed_wave' # Huge wave applied at the boundary
#scenario = 'slide'       # Slide wave form applied inside the domain

#------------------------------------------------------------------------------
# Filenames
#------------------------------------------------------------------------------
name_stem = cairns_dir+'cairns'
meshname =  name_stem + '.msh'

# Filename for locations where timeseries are to be produced
gauge_filename = cairns_dir + 'gauges.csv'

#------------------------------------------------------------------------------
# Domain definitions
#------------------------------------------------------------------------------
# bounding polygon for study area
bounding_polygon = anuga.read_polygon(cairns_dir+'extent.csv')

A = anuga.polygon_area(bounding_polygon) / 1000000.0
if __name__ == '_main_':
    print 'Area of bounding polygon = %.2f km^2' % A

#------------------------------------------------------------------------------
# Interior region definitions
#------------------------------------------------------------------------------
# Read interior polygons
poly_cairns = anuga.read_polygon(cairns_dir+'cairns.csv')
poly_island0 = anuga.read_polygon(cairns_dir+'islands.csv')
poly_island1 = anuga.read_polygon(cairns_dir+'islands1.csv')
poly_island2 = anuga.read_polygon(cairns_dir+'islands2.csv')
poly_island3 = anuga.read_polygon(cairns_dir+'islands3.csv')
poly_shallow = anuga.read_polygon(cairns_dir+'shallow.csv')

# Optionally plot points making up these polygons
#plot_polygons([bounding_polygon, poly_cairns, poly_island0, poly_island1,
#               poly_island2, poly_island3, poly_shallow],
#               style='boundingpoly', verbose=False)

# Define resolutions (max area per triangle) for each polygon
# Make these numbers larger to reduce the number of triangles in the model,
# and hence speed up the simulation

# bigger base_scale == less triangles
base_scale = 100000
base_scale = 1000000
default_res = 100 * base_scale   # Background resolution
islands_res = base_scale
cairns_res = base_scale
shallow_res = 5 * base_scale

# Define list of interior regions with associated resolutions
interior_regions = [[poly_cairns,  cairns_res],
                    [poly_island0, islands_res],
                    [poly_island1, islands_res],
                    [poly_island2, islands_res],
                    [poly_island3, islands_res],
                    [poly_shallow, shallow_res]]

#------------------------------------------------------------------------------
# Data for exporting ascii grid
#------------------------------------------------------------------------------
eastingmin = 363000
eastingmax = 418000
northingmin = 8026600
northingmax = 8145700

#------------------------------------------------------------------------------
# Data for landslide
#------------------------------------------------------------------------------
slide_origin = [451871, 8128376]   # Assume to be on continental shelf
slide_depth = 500.
