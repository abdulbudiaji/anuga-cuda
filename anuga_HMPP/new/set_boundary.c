#include "hmpp_fun.h"



#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp setBoundaryE codelet, target=CUDA args[*].transfer=atcall
#endif
void set_boundary_values_from_edges(
        int Nb,
        int N3,
        long vol_id[Nb],
        long edge_id[Nb],
        double boundary_values[Nb],
        double edge_values[N3])
{
    int k, id;
    
    /* Original function in quantity.py
        
        def set_boundary_values_from_edges(self):
            vol_id = self.domain.boundary_cells
            edge_id = self.domain.boundary_edges

            self.boundary_values =(self.edge_values.flat)[3*vol_id +edge_id]


        # Equals to below:


        for j in range(self.boundary_length):
            vol_id = self.domain.boundary_cells[j]
            edge_id = self.domain.boundary_edges[j]
            self.boundary_values[j] = self.edge_values[vol_id, edge_id]
    */

    for(k=0; k<Nb; k++)
    {
        id= 3*vol_id[k] + edge_id[k];

        boundary_values[k] = edge_values[id];
    }
}
