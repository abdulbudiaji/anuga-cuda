
// OpenHMPP ANUGA evolve function
//
// Zhe Weng 2013
#include <hmpp_fun.h>



int apply_fractional_steps( struct domain D)
{   return 0;}
int apply_protection_against_isolated_degenerate_timesteps(struct domain D)
{   return 0;}
int update_ghosts(struct domain D)
{   return 0;}
int update_extrema(struct domain D)
{   return 0;}
int log_operator_timestepping_statistics(struct domain D)
{   return 0;}



int update_boundary(struct domain D)
{   return 0;}



int update_timestep(struct domain D, double yieldstep, double finaltime)
{
    double timestep;

    apply_protection_against_isolated_degenerate_timesteps(D);
    
    timestep = fmin( D.CFL * D.flux_timestep, D.evolve_max_timestep);

    D.recorded_max_timestep = fmax( timestep, D.recorded_max_timestep);
    D.recorded_min_timestep = fmin( timestep, D.recorded_min_timestep);

    if ( timestep < D.evolve_min_timestep )
    {
        D.smallsteps += 1;
        
        if ( D.smallsteps > D.max_smallsteps)
        {
            D.smallsteps = 0;

            if ( D._order_ == 1)
            {
                timestep = D.evolve_min_timestep;
            } else
                D._order_ = 1;
        }
    } else {
        D.smallsteps = 0;
        if ( D._order_ == 1 && D.default_order == 2)
            D._order_ = 2;
    }

    if ( finaltime && D.time + timestep > finaltime)
        timestep = finaltime + D.time;

    if ( D.time + timestep > D.yieldtime)
        timestep = D.yieldtime - D.time;

    D.timestep = timestep;
    return 0;
}



int cmp(const void *a, const void *b)
{   
    return  *(double *)a - *(double *)b; 
}



double compute_fluxes(struct domain D)
{
    int i;
    switch ( D.compute_fluxes_method )
    {
        case 0: // 'original'
            break;
        case 1: // 'wb_1'
            break;
        case 2: // 'wb_2'
            #pragma hmpp cf_central callsite
            compute_fluxes_central_structure_CUDA(
                    D.number_of_elements, 
                    D.number_of_elements*3,
                    D.number_of_elements*6,
                    D.number_of_boundary_elements,

                    D.timestep_array,
                    D.neighbours,
                    D.neighbour_edges,
                    D.normals,
                    D.edgelengths,
                    D.radii,
                    D.areas,
                    D.tri_full_flag,
                    
                    D.stage_edge_values,
                    D.xmom_edge_values,
                    D.ymom_edge_values,
                    D.bed_edge_values,
                    
                    D.stage_boundary_values,
                    D.xmom_boundary_values,
                    D.ymom_boundary_values,
                    
                    D.stage_explicit_update,
                    D.xmom_explicit_update,
                    D.ymom_explicit_update,

                    D.max_speed,

                    D.evolve_max_timestep, 
                    D.g, 
                    D.epsilon,
                    D.h0,
                    D.limiting_threshold,
                    D.optimise_dry_cells);

            #pragma hmpp gravity callsite
            gravity_wb(
                D.number_of_elements,
                D.number_of_elements * 3,
                D.number_of_elements * 6,
                D.xmom_explicit_update,
                D.ymom_explicit_update,

                D.stage_vertex_values,
                D.stage_edge_values,
                D.stage_centroid_values,

                D.bed_edge_values,
                D.bed_centroid_values,

                D.vertex_coordinates,

                D.normals,
                D.areas,
                D.edgelengths,

                D.g
                );
            break;
        case 3: // 'wb_3'
            break;
        case 4: // 'tsunami'
            break;
        default:
            printf("unknown compute_fluxes_method\n");
            exit(EXIT_FAILURE);
            break;
    }

    // FIXME: needs transmit back
    // min
    //qsort(D.timestep, D.number_of_elements, sizeof(double), cmp);
    for (i=0;i < D.number_of_elements; i++)
    {
        if ( D.flux_timestep > D.timestep_array[i] )
            D.flux_timestep = D.timestep_array[i];
    }

    return D.flux_timestep;
}



int manning_friction_implicit(struct domain D)
{
    if ( D.use_sloped_mannings )
    {
        manning_friction_sloped(
                D.number_of_elements,
                D.g,
                D.minimum_allowed_height,

                D.vertex_coordinates,
                
                D.stage_centroid_values,
                D.bed_vertex_values,
                D.xmom_centroid_values,
                D.ymom_centroid_values,
                D.friction_centroid_values,

                D.xmom_semi_implicit_update,
                D.ymom_semi_implicit_update
                );
    } else {
        manning_friction_flat(
                D.number_of_elements,
                D.g,
                D.minimum_allowed_height,

                D.stage_centroid_values,
                D.bed_vertex_values,
                D.xmom_centroid_values,
                D.ymom_centroid_values,
                D.friction_centroid_values,

                D.xmom_semi_implicit_update,
                D.ymom_semi_implicit_update
                );
    }
    return 0;
}



int compute_forcing_terms(struct domain D)
{
    // FIXME:
    manning_friction_implicit(D);
    return 0;
}



int extrapolate_second_order_sw(struct domain D)
{
    if ( D.extrapolate_velocity_second_order )
    {
        // FIXME: OpenHMPP directives
        extrapolate_second_order_velocity_true(
                D.number_of_elements,
                D.minimum_allowed_height,

                D.stage_centroid_values,
                D.bed_centroid_values,
                D.xmom_centroid_values,
                D.xmom_centroid_store,
                D.ymom_centroid_values,
                D.ymom_centroid_store
                );

        // FIXME: OpenHMPP directives
        extrapolate_second_order_sw_true(
                D.number_of_elements,
                D.epsilon,
                D.minimum_allowed_height,
                D.beta_w,
                D.beta_w_dry,
                D.beta_uh,
                D.beta_uh_dry,
                D.beta_vh,
                D.beta_vh_dry,
                D.optimise_dry_cells,

                D.surrogate_neighbours,
                D.number_of_boundaries,
                
                D.centroid_coordinates,
                
                D.stage_centroid_values,
                D.bed_centroid_values,
                D.xmom_centroid_store,
                D.ymom_centroid_store,
                
                D.vertex_coordinates,

                D.stage_vertex_values,
                D.bed_vertex_values,
                D.xmom_vertex_values,
                D.ymom_vertex_values
                );
    } else {
        // FIXME: OpenHMPP directives
        extrapolate_second_order_sw_false(
                D.number_of_elements,
                D.epsilon,
                D.minimum_allowed_height,
                D.beta_w,
                D.beta_w_dry,
                D.beta_uh,
                D.beta_uh_dry,
                D.beta_vh,
                D.beta_vh_dry,
                D.optimise_dry_cells,

                D.surrogate_neighbours,
                D.number_of_boundaries,
                
                D.centroid_coordinates,
                
                D.stage_centroid_values,
                D.bed_centroid_values,
                D.xmom_centroid_values,
                D.ymom_centroid_values,
                
                D.vertex_coordinates,

                D.stage_vertex_values,
                D.bed_vertex_values,
                D.xmom_vertex_values,
                D.ymom_vertex_values
                );
    }
    return 0;
}



int update_conserved_quantities(struct domain D)
{
    //for name in self.conserved_quantities
    // stage
    update(
            D.number_of_elements,
            D.timestep,
            D.stage_centroid_values,
            D.stage_explicit_update,
            D.stage_semi_implicit_update
            );
    // FIXME: on device
    memset( D.stage_semi_implicit_update, 0, D.number_of_elements);
    // xmomentum
    update(
            D.number_of_elements,
            D.timestep,
            D.xmom_centroid_values,
            D.xmom_explicit_update,
            D.xmom_semi_implicit_update
            );
    // FIXME: on device
    memset( D.xmom_semi_implicit_update, 0, D.number_of_elements);
    // ymomentum
    update(
            D.number_of_elements,
            D.timestep,
            D.ymom_centroid_values,
            D.ymom_explicit_update,
            D.ymom_semi_implicit_update
            );
    // FIXME: on device
    memset( D.ymom_semi_implicit_update, 0, D.number_of_elements);
    return 0;
}



int backup_conserved_quantities(struct domain D)
{
    // FIXME: on device
    //for name in self.conserved_quantities
    // stage
    memcpy(D.stage_centroid_backup, D.stage_centroid_values, D.number_of_elements);
    // xmomentum
    memcpy(D.xmom_centroid_backup, D.xmom_centroid_values, D.number_of_elements);
    // ymomentum  centroid_
    memcpy(D.ymom_centroid_backup, D.ymom_centroid_values, D.number_of_elements);
    return 0;
}



int saxpy_conserved_quantities(struct domain D, double a, double b)
{
    //for name in self.conserved_quantities
    // stage
    saxpy_centroid_values(
            D.number_of_elements,
            a,
            b,
            D.stage_centroid_values,
            D.stage_centroid_backup
            );
    // xmomentum
    saxpy_centroid_values(
            D.number_of_elements,
            a,
            b,
            D.xmom_centroid_values,
            D.xmom_centroid_backup
            );
    // ymomentum
    saxpy_centroid_values(
            D.number_of_elements,
            a,
            b,
            D.ymom_centroid_values,
            D.ymom_centroid_backup
            );
    return 0;
}



int update_centroids_of_velocities_and_height( struct domain D)
{
    // elevation
    set_boundary_values_from_edges(
            D.number_of_boundary_elements,
            D.boundary_cells,
            D.boundary_edges,
            D.bed_boundary_values,
            D.bed_edge_values
            );

    _update_centroids_of_velocities_and_height(
            D.number_of_elements,
            D.number_of_boundary_elements,
            
            D.stage_centroid_values,
            D.xmom_centroid_values,
            D.ymom_centroid_values,
            D.height_centroid_values,
            D.bed_centroid_values,
            D.xvelocity_centroid_values,
            D.yvelocity_centroid_values,

            D.stage_boundary_values,
            D.xmom_boundary_values,
            D.ymom_boundary_values,
            D.height_boundary_values,
            D.bed_boundary_values,
            D.xvelocity_boundary_values,
            D.yvelocity_boundary_values
            );
    return 0;
}



int update_other_quantities( struct domain D )
{
    // 'yusuke'
    if ( D.flow_algorithm == 2)
        return 0;

    update_centroids_of_velocities_and_height(D);

    //for name in ['height', 'xvelocity', 'yvelocity']
    // height
    extrapolate_first_order(
            D.number_of_elements,
            D.height_centroid_values,
            D.height_edge_values,
            D.height_vertex_values
            );
    // xvelocity
    extrapolate_first_order(
            D.number_of_elements,
            D.xvelocity_centroid_values,
            D.xvelocity_edge_values,
            D.xvelocity_vertex_values
            );
    // yvelocity
    extrapolate_first_order(
            D.number_of_elements,
            D.yvelocity_centroid_values,
            D.yvelocity_edge_values,
            D.yvelocity_vertex_values
            );
    return 0;
}



int protect_against_infinitesimal_and_negative_heights( struct domain D)
{
    // 'tsunami'
    if ( D.flow_algorithm == 1)
    {
        protect_swb2(
            D.number_of_elements,
            D.number_of_elements * 3,
            D.minimum_allowed_height,
            D.maximum_allowed_speed,
            D.epsilon,
            
            D.stage_centroid_values,
            D.stage_vertex_values,
            
            D.bed_centroid_values,
            D.bed_vertex_values,

            D.xmom_centroid_values,
            D.ymom_centroid_values,
            D.areas
            );
    } else {
        protect_sw(
            D.number_of_elements,
            D.number_of_elements * 3,
            D.minimum_allowed_height,
            D.maximum_allowed_speed,
            D.epsilon,
            
            D.stage_centroid_values,
            D.bed_centroid_values,
            D.xmom_centroid_values,
            D.ymom_centroid_values
            );
    }
    return 0;
}



int distribute_using_edge_limiter(struct domain D)
{
    return 0;
}

int distribute_using_vertex_limiter(struct domain D)
{
    protect_against_infinitesimal_and_negative_heights(D);

    if ( D.optimised_gradient_limiter )
    {
        if ( D._order_ == 1)
        {
            // for name in self.conserved_quantities:
            // stage
            extrapolate_first_order(
                    D.number_of_elements,
                    D.stage_centroid_values,
                    D.stage_edge_values,
                    D.stage_vertex_values
                    );
            // xmomentum
            extrapolate_first_order(
                    D.number_of_elements,
                    D.xmom_centroid_values,
                    D.xmom_edge_values,
                    D.xmom_vertex_values
                    );
            // ymomentum
            extrapolate_first_order(
                    D.number_of_elements,
                    D.ymom_centroid_values,
                    D.ymom_edge_values,
                    D.ymom_vertex_values
                    );
        }
        else if (D._order_ == 2)
        {
            extrapolate_second_order_sw( D );
        }
        else
        {
            printf("Error: unknown order\n");
            exit(EXIT_FAILURE);
        }
    } else {
        // for name in self.conserved_quantities:
        if ( D._order_ == 1)
        {
            // for name in self.conserved_quantities:
            // stage
            extrapolate_first_order(
                    D.number_of_elements,
                    D.stage_centroid_values,
                    D.stage_edge_values,
                    D.stage_vertex_values
                    );
            // xmomentum
            extrapolate_first_order(
                    D.number_of_elements,
                    D.xmom_centroid_values,
                    D.xmom_edge_values,
                    D.xmom_vertex_values
                    );
            // ymomentum
            extrapolate_first_order(
                    D.number_of_elements,
                    D.ymom_centroid_values,
                    D.ymom_edge_values,
                    D.ymom_vertex_values
                    );
        }
        else if (D._order_ == 2)
        {
            // stage
            extrapolate_second_order_and_limit_by_vertex(
                    D.number_of_elements,
                    D.stage_beta,

                    D.centroid_coordinates,
                    D.vertex_coordinates,

                    D.number_of_boundaries,
                    D.surrogate_neighbours,
                    D.neighbours,

                    D.stage_centroid_values,
                    D.stage_vertex_values,
                    D.stage_edge_values,
                    D.stage_x_gradient,
                    D.stage_y_gradient
                    );
            // xmomentum
            extrapolate_second_order_and_limit_by_vertex(
                    D.number_of_elements,
                    D.xmom_beta,

                    D.centroid_coordinates,
                    D.vertex_coordinates,

                    D.number_of_boundaries,
                    D.surrogate_neighbours,
                    D.neighbours,

                    D.xmom_centroid_values,
                    D.xmom_vertex_values,
                    D.xmom_edge_values,
                    D.xmom_x_gradient,
                    D.xmom_y_gradient
                    );
            // ymomentum
            extrapolate_second_order_and_limit_by_vertex(
                    D.number_of_elements,
                    D.ymom_beta,

                    D.centroid_coordinates,
                    D.vertex_coordinates,

                    D.number_of_boundaries,
                    D.surrogate_neighbours,
                    D.neighbours,

                    D.ymom_centroid_values,
                    D.ymom_vertex_values,
                    D.ymom_edge_values,
                    D.ymom_x_gradient,
                    D.ymom_y_gradient
                    );
        }
        else
        {
            printf("Error: unknown order\n");
            exit(EXIT_FAILURE);
        }
        
    }

    #pragma hmpp balance callsite 
    balance_deep_and_shallow(
           D.number_of_elements,
           D.H0,
           D.alpha_balance,
           D.tight_slope_limiters,
           D.use_centroid_velocities,

           D.stage_centroid_values,
           D.bed_centroid_values,
           D.stage_vertex_values,
           D.bed_vertex_values,

           D.xmom_centroid_values,
           D.ymom_centroid_values,

           D.xmom_vertex_values,
           D.ymom_vertex_values
           );
    // for name in self.conserved_quantities:
    // stage
    interpolate_from_vertices_to_edges(
            D.number_of_elements,
            D.stage_vertex_values,
            D.stage_edge_values
            );
    // xmomentum
    interpolate_from_vertices_to_edges(
            D.number_of_elements,
            D.xmom_vertex_values,
            D.xmom_edge_values
            );
    // ymomentum
    interpolate_from_vertices_to_edges(
            D.number_of_elements,
            D.ymom_vertex_values,
            D.ymom_edge_values
            );
    return 0;
}



int distribute_to_vertices_and_edges(struct domain D)
{
    // FIXME: compute_fluxes_method == 'tsunami'
    if ( D.compute_fluxes_method == 1)
    {
        // FIXME : OpenHMPP directives
        protect_swb2(
            D.number_of_elements,
            D.number_of_elements * 3,
            D.minimum_allowed_height,
            D.maximum_allowed_speed,
            D.epsilon,
            
            D.stage_centroid_values,
            D.stage_vertex_values,

            D.bed_centroid_values,
            D.bed_vertex_values,
            
            D.xmom_centroid_values,
            D.ymom_centroid_values,
            D.areas
            );


        // FIXME : OpenHMPP directives
        // from swb2_domain_ext
        extrapolate_second_order_edge_sw(
            D.number_of_elements,
            D.optimise_dry_cells,
            D.extrapolate_velocity_second_order,
            
            D.epsilon,
            D.minimum_allowed_height,
            D.beta_w,
            D.beta_w_dry,
            D.beta_uh,
            D.beta_uh_dry,
            D.beta_vh,
            D.beta_vh_dry,

            D.surrogate_neighbours,
            D.number_of_boundaries,
            
            D.centroid_coordinates,
            
            D.stage_centroid_values,
            D.xmom_centroid_values,
            D.ymom_centroid_values,
            D.bed_centroid_values,
            
            D.edge_coordinates,
            
            D.stage_edge_values,
            D.xmom_edge_values,
            D.ymom_edge_values,
            D.bed_edge_values,

            D.stage_vertex_values,
            D.xmom_vertex_values,
            D.ymom_vertex_values,
            D.bed_vertex_values,

            D.stage_centroid_store,
            D.xmom_centroid_store,
            D.ymom_centroid_store,

            D.min_bed_edge_values,
            D.max_bed_edge_values,
            D.count_wet_neighbours
            );
    }
    else if ( D.use_edge_limiter )
        distribute_using_edge_limiter(D);
    else
        distribute_using_vertex_limiter(D);
    return 0;
}



int evolve_one_euler_step(struct domain D, double yieldstep, double finaltime)
{
    compute_fluxes( D );
    compute_forcing_terms( D );
    update_timestep( D, yieldstep, finaltime);
    update_conserved_quantities( D );
    return 0;
}



int evolve_one_rk2_step(struct domain D, double yieldstep, double finaltime)
{
    return 0;
}
int evolve_one_rk3_step(struct domain D, double yieldstep, double finaltime)
{
    return 0;
}


int evolve( struct domain D, 
            double yieldstep, 
            double finaltime,
            double duration,
            double epsilon,
            int skip_initial_step
            )
{
    double initial_time;


    assert( D.beta_w >= 0 && D.beta_w <= 2.0 );


    distribute_to_vertices_and_edges(D);

    if (D.time != D.starttime)
        D.time = D.starttime;
    if ( ! yieldstep)
        yieldstep = D.evolve_max_timestep;

    D._order_ = D.default_order;

    assert( finaltime >= D.starttime );

    if (finaltime)
        D.finaltime = finaltime;
    else if (duration)
        D.finaltime = duration;
    else
    {
        printf("Only one of finaltime and duration may be specified\n");
        exit(EXIT_FAILURE);
    }

    D.yieldtime = D.time + yieldstep;

    D.recorded_min_timestep = D.evolve_max_timestep;
    D.recorded_max_timestep = D.evolve_min_timestep;
    D.number_of_steps = 0;
    D.number_of_first_order_steps = 0;

    update_ghosts(D);
    distribute_to_vertices_and_edges(D);
    update_boundary(D);

    update_extrema(D);        

    if ( !skip_initial_step )
        return D.time;

    while(1)
    {
        initial_time = D.time;

        // euler
        if ( D.timestepping_method == 1)
            evolve_one_euler_step(D, yieldstep, D.finaltime);
        // rk2
        else if (D.timestepping_method == 2)
            evolve_one_rk2_step(D, yieldstep, D.finaltime);
        // rk3
        else if (D.timestepping_method == 3)
            evolve_one_rk3_step(D, yieldstep, D.finaltime);

        apply_fractional_steps(D);

        D.time = initial_time + D.timestep;

        update_ghosts(D);

        distribute_to_vertices_and_edges(D);

        update_boundary(D);

        update_other_quantities(D);

        update_extrema(D);

        D.number_of_steps += 1;
        if ( D._order_ == 1 )
            D.number_of_first_order_steps += 1;

        if (D.finaltime && D.time >= D.finaltime - epsilon)
        {
            if (D.time > D.finaltime)
            {    
                printf("WARNING (domain.py): time overshot finaltime.\n");
                exit(EXIT_FAILURE);
            }
            D.time = D.finaltime;
            log_operator_timestepping_statistics(D);
            
            // FIXME: exit the whole program
            return D.time;
        }

        if (D.time >= D.yieldtime)
        {
            log_operator_timestepping_statistics(D);
            // FIXME: yield(self.get_time())
            return D.time;

            D.yieldtime += yieldstep;
            D.recorded_min_timestep = D.evolve_max_timestep;
            D.recorded_max_timestep = D.evolve_min_timestep;
            D.number_of_steps = 0;
            D.number_of_first_order_steps = 0;
            memset(D.max_speed, 0, D.number_of_elements);
        }
    }



}
