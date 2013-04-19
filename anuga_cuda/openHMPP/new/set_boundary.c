void set_boundary_values_from_edges(
        int N,
        long * vol_id,
        long * edge_id,
        double * boundary_values,
        double * edge_values)
{
    int k, id;
    
    for(k=0; k<N; k++)
    {
        id= 3*vol_id[k] + edge_id[k];

        boundary_values[k] = edge_values[id];
    }
}
