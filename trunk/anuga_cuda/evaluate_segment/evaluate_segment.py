
import numpy
from pycuda import driver as drv
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import \
        Dirichlet_boundary
from anuga.shallow_water.boundaries import Reflective_boundary

from anuga_cuda import generate_merimbula_domain
from anuga_cuda import generate_cairns_domain

domain1 = generate_merimbula_domain( gpu=False )
domain2 = generate_merimbula_domain( gpu=True )
domain2.equip_kernel_functions()


domain1.distribute_to_vertices_and_edges()
domain2.distribute_to_vertices_and_edges()
        

#domain1.update_boundary()
for tag in domain1.tag_boundary_cells:
    B1 = domain1.boundary_map[tag]
    B2 = domain2.boundary_map[tag]
    if B1 is None:
        continue
    else:
        B1.evaluate_segment(domain1, domain1.tag_boundary_cells[tag] )


        ids = domain2.tag_boundary_cells[tag]
        N = len(ids)
        W1 = 32
        W2 = 1
        W3 = 1
        
        if isinstance( B2, Reflective_boundary):
            domain2.evaluate_segment_reflective_func(
                numpy.int32( N ),
                drv.In( numpy.asarray(ids) ),
                drv.In( numpy.asarray( domain2.boundary_cells[ids]) ),
                drv.In( numpy.asarray( domain2.boundary_edges[ids]) ),

                drv.In( domain2.normals ),
                drv.In( domain2.quantities['stage'].edge_values),
                drv.In( domain2.quantities['elevation'].edge_values),
                drv.In( domain2.quantities['height'].edge_values),
                drv.In( domain2.quantities['xmomentum'].edge_values),
                drv.In( domain2.quantities['ymomentum'].edge_values),
                drv.In( domain2.quantities['xvelocity'].edge_values),
                drv.In( domain2.quantities['yvelocity'].edge_values),
        
                drv.InOut( domain2.quantities['stage'].boundary_values),
                drv.InOut( domain2.quantities['elevation'].boundary_values),
                drv.InOut( domain2.quantities['height'].boundary_values),
                drv.InOut( domain2.quantities['xmomentum'].boundary_values),
                drv.InOut( domain2.quantities['ymomentum'].boundary_values),
                drv.InOut( domain2.quantities['xvelocity'].boundary_values),
                drv.InOut( domain2.quantities['yvelocity'].boundary_values),

                block = (W1, W2, W3),
                grid = ( (N+ W1*W2*W3 - 1)/(W1*W2*W3), 1)
                )

            s1 = domain1.quantities['stage'].boundary_values
            e1 = domain1.quantities['elevation'].boundary_values
            h1 = domain1.quantities['height'].boundary_values
            xm1 = domain1.quantities['xmomentum'].boundary_values
            ym1 = domain1.quantities['ymomentum'].boundary_values
            xv1 = domain1.quantities['xvelocity'].boundary_values
            yv1 = domain1.quantities['yvelocity'].boundary_values


            s2 = domain2.quantities['stage'].boundary_values
            e2 = domain2.quantities['elevation'].boundary_values
            h2 = domain2.quantities['height'].boundary_values
            xm2 = domain2.quantities['xmomentum'].boundary_values
            ym2 = domain2.quantities['ymomentum'].boundary_values
            xv2 = domain2.quantities['xvelocity'].boundary_values
            yv2 = domain2.quantities['yvelocity'].boundary_values


            print  "%s %d -- Reflective_boundary" % (tag, N)
            cnt_s = 0
            cnt_e = 0
            cnt_h = 0
            cnt_xm = 0
            cnt_ym = 0
            cnt_xv = 0
            cnt_yv = 0
            for i in range(N):
                k = ids[i]
                if s1[k] != s2[k]:
                    cnt_s += 1
                if e1[k] != e2[k]:
                    cnt_e += 1
                if h1[k] != h2[k]:
                    cnt_h += 1
                if xm1[k] != xm2[k]:
                    cnt_xm += 1
                if ym1[k] != ym2[k]:
                    cnt_ym += 1
                if xv1[k] != xv2[k]:
                    cnt_xv += 1
                if yv1[k] != yv2[k]:
                    cnt_yv += 1
            print "    ",cnt_s, cnt_e, cnt_h, cnt_xm, cnt_ym, cnt_xv, cnt_yv


        elif isinstance( B2, Dirichlet_boundary):
            
            print "%s  %d -- Dirichlet_boundary" % (tag, N)
            
            q_bdry = B2.dirichlet_values
            conserved_quantities = True
            if len(q_bdry) == len(domain2.evolved_quantities):
                conserved_quantities = False
            if conserved_quantities:
                for j, name in enumerate( domain2.evolved_quantities):
                    Q2 = domain2.quantities[name]
                    self.evaluate_segment_dirichlet_1_func(
                        numpy.int32(N),
                        drv.In( numpy.asarray(ids) ),
                        drv.In( 
                            numpy.asarray( domain2.boundary_cells[ids]) ),
                        drv.In( 
                            numpy.asarray( domain2.boundary_edges[ids]) ),

                        drv.InOut( Q2.boundary_values ),
                        drv.In( Q2.edge_values),
                        block = (W1, W2, W3),
                        grid = ( (N+ W1*W2*W3 - 1)/(W1*W2*W3), 1)
                        )

                    cnt = 0
                    Q1 = domain1.quantities[name].boundary_values
                    for i in range(N):
                        k = ids[i]
                        if Q1[k] != Q2[k]:
                            cnt += 1
                            if cnt < 5 :
                                print "        %d, %f, %f"%(k, Q1[k], Q2[k])
                    print "    %s, %d" % (name, cnt)

            if conserved_quantities:
                quantities = domain2.conserved_quantities
            else:
                quantities = domain2.evolved_quantities



            for j, name in enumerate(quantities):
                Q2 = domain2.quantities[name].boundary_values
                domain2.evaluate_segment_dirichlet_2_func(
                    numpy.int32(N),
                    numpy.float64(q_bdry[j]),
                    drv.In( numpy.asanyarray( ids ) ),
                    drv.InOut( Q2 ),
                    block = (W1, W2, W3),
                    grid = ( (N+ W1*W2*W3 - 1)/(W1*W2*W3), 1)
                    )

                cnt = 0
                Q1 = domain1.quantities[name].boundary_values
                for i in range(N):
                    k = ids[i]
                    if Q1[k] != Q2[k]:
                        cnt += 1
                        if cnt < 5 :
                            print "        %d, %f, %f" % (k, Q1[k], Q2[k])
                print "    %s, %d" % (name, cnt)


        else:
            raise Exception("wrong type")
                    


