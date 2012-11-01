
// Python - C extension module for gpu_domain.py
//
// To compile (Python2.3):
//  gcc -c %.c -I/usr/include/python2.3 -o %.o -Wall -O
//  gcc -shared %.o  -o %.so cudafun.o -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lm
//
// or use python compile.py
//
// See the module shallow_water_domain.py for more documentation on 
// how to use this module
//
// Includes Cuda GPU usage: 
// memory allocation in first iteration, copying, and deletion in final
// iteration;
// requires linkage of cudafun.o, CUDA library
//


#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "numpy_shim.h"
#include "util_ext.h"

#include "cudafun.h"

//=========================================================================
// Python Glue
//=========================================================================

PyObject *compute_fluxes_ext_central_new_gpu(PyObject *self, PyObject *args) {
  /*Compute all fluxes and the timestep suitable for all volumes
    in domain. GPU version.

    Compute total flux for each conserved quantity using "flux_function_central"

    Fluxes across each edge are scaled by edgelengths and summed up
    Resulting flux is then scaled by area and stored in
    explicit_update for each of the three conserved quantities
    stage, xmomentum and ymomentum

    The maximal allowable speed computed by the flux_function for each volume
    is converted to a timestep that must not be exceeded. The minimum of
    those is computed as the next overall timestep.

	First set up data structures on GPU, then copy data, call computational kernel,
	copy back.

    Python call:
    timestep = compute_fluxes(timestep, domain, stage, xmom, ymom, bed)


    Post conditions:
      domain.explicit_update is reset to computed flux values
      returns timestep which is the largest step satisfying all volumes.


  */
    PyObject 
    *domain,
    *stage, 
    *xmom, 
    *ymom, 
    *bed;

    PyArrayObject 
    *neighbours, 
    *neighbour_edges,
    *normals, 
    *edgelengths, 
    *radii, 
    *areas,
    *tri_full_flag,
    *stage_edge_values,
    *xmom_edge_values,
    *ymom_edge_values,
    *bed_edge_values,
    *stage_boundary_values,
    *xmom_boundary_values,
    *ymom_boundary_values,
    *stage_explicit_update,
    *xmom_explicit_update,
    *ymom_explicit_update,
    *already_computed_flux, //Tracks whether the flux across an edge has already been computed
    *max_speed_array; //Keeps track of max speeds for each triangle

    static int iteration = 0;
    // TODO: get final_iteration
    int final_iter = 1;
    double timestep, epsilon, H0, g;
    int i,optimise_dry_cells = 0;
    
    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "dOOOOO", &timestep, &domain, &stage, &xmom, &ymom, &bed )) {
      report_python_error(AT, "could not parse input arguments");
      return NULL;
    }

    epsilon           = get_python_double(domain,"epsilon");
    H0                = get_python_double(domain,"H0");
    g                 = get_python_double(domain,"g");
    //optimise_dry_cells = get_python_integer(domain,"optimse_dry_cells");
    
    neighbours        = get_consecutive_array(domain, "neighbours");
    neighbour_edges   = get_consecutive_array(domain, "neighbour_edges"); 
    normals           = get_consecutive_array(domain, "normals");
    edgelengths       = get_consecutive_array(domain, "edgelengths");    
    radii             = get_consecutive_array(domain, "radii");    
    areas             = get_consecutive_array(domain, "areas");    
    tri_full_flag     = get_consecutive_array(domain, "tri_full_flag");
    already_computed_flux  = get_consecutive_array(domain, "already_computed_flux");
    max_speed_array   = get_consecutive_array(domain, "max_speed");
    
    stage_edge_values = get_consecutive_array(stage, "edge_values");    
    xmom_edge_values  = get_consecutive_array(xmom, "edge_values");    
    ymom_edge_values  = get_consecutive_array(ymom, "edge_values"); 
    bed_edge_values   = get_consecutive_array(bed, "edge_values");    

    stage_boundary_values = get_consecutive_array(stage, "boundary_values");    
    xmom_boundary_values  = get_consecutive_array(xmom, "boundary_values");    
    ymom_boundary_values  = get_consecutive_array(ymom, "boundary_values"); 

    stage_explicit_update = get_consecutive_array(stage, "explicit_update");    
    xmom_explicit_update  = get_consecutive_array(xmom, "explicit_update");    
    ymom_explicit_update  = get_consecutive_array(ymom, "explicit_update"); 

  int number_of_elements = stage_edge_values -> dimensions[0];

  printf("glue: dt: %f,  noe: %i, eps: %f, H0: %f, g: %f\n",
	timestep, number_of_elements, epsilon, H0, g); fflush(stdout);

  // get array dimensions
  size_t number_of_neighbours        = neighbours -> dimensions[0];
  size_t number_of_neighbour_edges   = neighbour_edges -> dimensions[0];
  size_t number_of_normals           = normals -> dimensions[0];
  size_t number_of_edgelengths       = edgelengths -> dimensions[0];
  size_t number_of_radii             = radii -> dimensions[0];
  size_t number_of_areas             = areas -> dimensions[0];
  size_t number_of_tri_full_flag     = tri_full_flag -> dimensions[0];
  //size_t number_of_already_computed  = already_computed_flux -> dimensions[0];
  size_t number_of_max_speed_array   = max_speed_array -> dimensions[0];
    
  size_t number_of_stage_edge_values = stage_edge_values -> dimensions[0];
  size_t number_of_xmom_edge_values  = xmom_edge_values -> dimensions[0];
  size_t number_of_ymom_edge_values  = ymom_edge_values -> dimensions[0];
  size_t number_of_bed_edge_values   = bed_edge_values -> dimensions[0];

  size_t number_of_stage_boundary_values = stage_boundary_values -> dimensions[0];
  size_t number_of_xmom_boundary_values  = xmom_boundary_values -> dimensions[0];
  size_t number_of_ymom_boundary_values  = ymom_boundary_values -> dimensions[0];
  
  size_t number_of_stage_explicit_update = stage_explicit_update -> dimensions[0];
  size_t number_of_xmom_explicit_update  = xmom_explicit_update -> dimensions[0];
  size_t number_of_ymom_explicit_update  = ymom_explicit_update -> dimensions[0];


   // extract C arrays from python wrapper
   long*    c_neighbours         	 =  (long*) neighbours -> data;
   long*    c_neighbour_edges    	 =  (long*) neighbour_edges -> data;
   double*  c_normals            	 =  (double*) normals -> data;
   double*  c_edgelengths        	 =  (double*) edgelengths -> data; 
   double*  c_radii              	 =  (double*) radii -> data; 
   double*  c_areas              	 =  (double*) areas -> data;
   long*    c_tri_full_flag       	 =  (long*) tri_full_flag -> data;
   //long*    c_already_computed_flux	 =  (long*) already_computed_flux -> data;
   double*  c_max_speed_array    	 =  (double*) max_speed_array -> data;

   double*  c_stage_edge_values  =  (double*) stage_edge_values -> data;
   double*  c_xmom_edge_values   =  (double*) xmom_edge_values -> data;
   double*  c_ymom_edge_values   =  (double*) ymom_edge_values -> data;
   double*  c_bed_edge_values    =  (double*) bed_edge_values -> data;
            
   double*  c_stage_boundary_values  =  (double*) stage_boundary_values -> data;
   double*  c_xmom_boundary_values   =  (double*) xmom_boundary_values -> data;
   double*  c_ymom_boundary_values   =  (double*) ymom_boundary_values -> data;
            
   double*  c_stage_explicit_update  =  (double*) stage_explicit_update -> data;
   double*  c_xmom_explicit_update   =  (double*) xmom_explicit_update -> data;
   double*  c_ymom_explicit_update   =  (double*) ymom_explicit_update -> data;

	setKernelDims( 32, 96 );
	printKernelDims();

   static long*    gpu_neighbours         	  = NULL; 
   static long*    gpu_neighbour_edges    	  = NULL; 
   static double*  gpu_normals            	  = NULL;
   static double*  gpu_edgelengths        	  = NULL;
   static double*  gpu_radii              	  = NULL;
   static double*  gpu_areas              	  = NULL;
   static long*    gpu_tri_full_flag           = NULL;
   //static long*    gpu_already_computed_flux  = NULL;
   static double*  gpu_max_speed_array    	  = NULL;
                                                                                         
   static double*  gpu_stage_edge_values  = NULL;
   static double*  gpu_xmom_edge_values   = NULL;
   static double*  gpu_ymom_edge_values   = NULL;
   static double*  gpu_bed_edge_values    = NULL;
                                                                                              
   static double*  gpu_stage_boundary_values = NULL;
   static double*  gpu_xmom_boundary_values  = NULL;
   static double*  gpu_ymom_boundary_values  = NULL;
                                                                                                      
   static double*  gpu_stage_explicit_update = NULL;
   static double*  gpu_xmom_explicit_update  = NULL;
   static double*  gpu_ymom_explicit_update  = NULL;
	
   
	if( 0 == iteration ) {
	selectDevice(0,1,"");
   // allocate GPU arrays
   gpu_neighbours         	  = (long*)  	allocDeviceMemory( number_of_neighbours       * sizeof(long) ); 
   gpu_neighbour_edges    	  = (long*)  	allocDeviceMemory( number_of_neighbour_edges  * sizeof(long) ); 
   gpu_normals            	  = (double*)	allocDeviceMemory( number_of_normals          * sizeof(double) );
   gpu_edgelengths        	  = (double*)	allocDeviceMemory( number_of_edgelengths      * sizeof(double) );
   gpu_radii              	  = (double*)	allocDeviceMemory( number_of_radii            * sizeof(double) );
   gpu_areas              	  = (double*)	allocDeviceMemory( number_of_areas            * sizeof(double) );
   gpu_tri_full_flag           = (long*)  	allocDeviceMemory( number_of_tri_full_flag    * sizeof(long) );
   //gpu_already_computed_flux = (long*)  	allocDeviceMemory( number_of_already_computed * sizeof(long) );
   gpu_max_speed_array    	  = (double*)	allocDeviceMemory( number_of_max_speed_array  * sizeof(double) );
                                                                         
   gpu_stage_edge_values  = (double*) allocDeviceMemory(    number_of_stage_edge_values   * sizeof(double) );
   gpu_xmom_edge_values   = (double*) allocDeviceMemory(    number_of_xmom_edge_values    * sizeof(double) );
   gpu_ymom_edge_values   = (double*) allocDeviceMemory(    number_of_ymom_edge_values    * sizeof(double) );
   gpu_bed_edge_values    = (double*) allocDeviceMemory(    number_of_bed_edge_values     * sizeof(double) );
                                                                              
   gpu_stage_boundary_values = (double*) allocDeviceMemory( number_of_stage_boundary_values * sizeof(double) );
   gpu_xmom_boundary_values  = (double*) allocDeviceMemory( number_of_xmom_boundary_values  * sizeof(double) );
   gpu_ymom_boundary_values  = (double*) allocDeviceMemory( number_of_ymom_boundary_values  * sizeof(double) );
                                                                                      
   gpu_stage_explicit_update = (double*) allocDeviceMemory(  number_of_stage_explicit_update* sizeof(double) );
   gpu_xmom_explicit_update  = (double*) allocDeviceMemory(  number_of_xmom_explicit_update  * sizeof(double) );
   gpu_ymom_explicit_update  = (double*) allocDeviceMemory(  number_of_ymom_explicit_update  * sizeof(double) );

   // Constant quantities copied only in first iteration	
   copyHostToDevice( gpu_neighbours         	 , c_neighbours         	 , number_of_neighbours       * sizeof(long) ); 
   copyHostToDevice( gpu_neighbour_edges    	 , c_neighbour_edges    	 , number_of_neighbour_edges  * sizeof(long) ); 
   copyHostToDevice( gpu_normals            	 , c_normals            	 , number_of_normals          * sizeof(double) );
   copyHostToDevice( gpu_edgelengths        	 , c_edgelengths        	 , number_of_edgelengths      * sizeof(double) );
   copyHostToDevice( gpu_radii              	 , c_radii              	 , number_of_radii            * sizeof(double) );
   copyHostToDevice( gpu_areas              	 , c_areas              	 , number_of_areas            * sizeof(double) );
   }

   copyHostToDevice( gpu_tri_full_flag        ,    c_tri_full_flag        ,  number_of_tri_full_flag    * sizeof(long) );
   //copyHostToDevice( gpu_already_computed_flux,    c_already_computed_flux,  number_of_already_computed * sizeof(long) );
   //copyHostToDevice( gpu_max_speed_array    	 , c_max_speed_array    	 , number_of_max_speed_array  * sizeof(double) );
                                                                                                        
   copyHostToDevice( gpu_stage_edge_values  ,    c_stage_edge_values  ,    number_of_stage_edge_values   * sizeof(double) );
   copyHostToDevice( gpu_xmom_edge_values   ,    c_xmom_edge_values   ,    number_of_xmom_edge_values    * sizeof(double) );
   copyHostToDevice( gpu_ymom_edge_values   ,    c_ymom_edge_values   ,    number_of_ymom_edge_values    * sizeof(double) );
   copyHostToDevice( gpu_bed_edge_values    ,    c_bed_edge_values    ,    number_of_bed_edge_values     * sizeof(double) );
                                                                                                             
   copyHostToDevice( gpu_stage_boundary_values  ,c_stage_boundary_values  ,number_of_stage_boundary_values * sizeof(double) );
   copyHostToDevice( gpu_xmom_boundary_values   ,c_xmom_boundary_values   ,number_of_xmom_boundary_values  * sizeof(double) );
   copyHostToDevice( gpu_ymom_boundary_values   ,c_ymom_boundary_values   ,number_of_ymom_boundary_values  * sizeof(double) );
                                                                                                             
   /*copyHostToDevice( gpu_stage_explicit_update  ,c_stage_explicit_update  ,number_of_stage_explicit_update* sizeof(double) );
   copyHostToDevice( gpu_xmom_explicit_update   ,c_xmom_explicit_update   ,number_of_xmom_explicit_update  * sizeof(double) );
   copyHostToDevice( gpu_ymom_explicit_update   ,c_ymom_explicit_update   ,number_of_ymom_explicit_update  * sizeof(double) );*/
	
  // initialize explicit updates to zero (possibly superfluous)
  _set_to_default( gpu_stage_explicit_update, gpu_xmom_explicit_update, gpu_ymom_explicit_update, number_of_stage_explicit_update, 0.0 );

  
  // Call underlying flux computation routine and update 
  // the explicit update arrays 
  timestep = _compute_fluxes_central(number_of_elements,
                     timestep,
                     epsilon,
                     H0,
                     g,
                     gpu_neighbours         	 ,
                     gpu_neighbour_edges    	 ,
                     gpu_normals            	 ,
                     gpu_edgelengths        	 ,
                     gpu_radii              	 ,
                     gpu_areas              	 ,
                     gpu_tri_full_flag        ,  
				     gpu_stage_edge_values  ,    
                     gpu_xmom_edge_values   ,    
                     gpu_ymom_edge_values   ,    
                     gpu_bed_edge_values    ,    
                     gpu_stage_boundary_values  ,
                     gpu_xmom_boundary_values   ,
                     gpu_ymom_boundary_values   ,
                     gpu_stage_explicit_update  ,
                     gpu_xmom_explicit_update   ,
                     gpu_ymom_explicit_update   ,
                     //gpu_already_computed_flux  ,
                     gpu_max_speed_array        ,
					 optimise_dry_cells);
 
/*                     (long*) neighbours -> data,
                     (long*) neighbour_edges -> data,
                     (double*) normals -> data,
                     (double*) edgelengths -> data, 
                     (double*) radii -> data, 
                     (double*) areas -> data,
                     (long*) tri_full_flag -> data,
                     (double*) stage_edge_values -> data,
                     (double*) xmom_edge_values -> data,
                     (double*) ymom_edge_values -> data,
                     (double*) bed_edge_values -> data,
                     (double*) stage_boundary_values -> data,
                     (double*) xmom_boundary_values -> data,
                     (double*) ymom_boundary_values -> data,
                     (double*) stage_explicit_update -> data,
                     (double*) xmom_explicit_update -> data,
                     (double*) ymom_explicit_update -> data,
                     (long*) already_computed_flux -> data,
                     (double*) max_speed_array -> data,
                     optimise_dry_cells);
*/

   // copy GPU to Host memory	
   /*copyDeviceToHost( c_neighbours         	 , gpu_neighbours         	 , number_of_neighbours       * sizeof(long) ); 
   copyDeviceToHost( c_neighbour_edges    	 , gpu_neighbour_edges    	 , number_of_neighbour_edges  * sizeof(long) ); 
   copyDeviceToHost( c_normals            	 , gpu_normals            	 , number_of_normals          * sizeof(double) );
   copyDeviceToHost( c_edgelengths        	 , gpu_edgelengths        	 , number_of_edgelengths      * sizeof(double) );
   copyDeviceToHost( c_radii              	 , gpu_radii              	 , number_of_radii            * sizeof(double) );
   copyDeviceToHost( c_areas              	 , gpu_areas              	 , number_of_areas            * sizeof(double) );
   copyDeviceToHost( c_tri_full_flag        ,  gpu_tri_full_flag        ,      number_of_tri_full_flag    * sizeof(long) );
   //copyDeviceToHost( c_already_computed_flux,  gpu_already_computed_flux,      number_of_already_computed * sizeof(long) );
   copyDeviceToHost( c_max_speed_array    	 , gpu_max_speed_array    	 , number_of_max_speed_array  * sizeof(double) );
                                                                                                        
   copyDeviceToHost( c_stage_edge_values  ,    gpu_stage_edge_values  ,        number_of_stage_edge_values   * sizeof(double) );
   copyDeviceToHost( c_xmom_edge_values   ,    gpu_xmom_edge_values   ,        number_of_xmom_edge_values    * sizeof(double) );
   copyDeviceToHost( c_ymom_edge_values   ,    gpu_ymom_edge_values   ,        number_of_ymom_edge_values    * sizeof(double) );
   copyDeviceToHost( c_bed_edge_values    ,    gpu_bed_edge_values    ,        number_of_bed_edge_values     * sizeof(double) );
                                                                                                             
   copyDeviceToHost( c_stage_boundary_values  ,gpu_stage_boundary_values  ,    number_of_stage_boundary_values * sizeof(double) );
   copyDeviceToHost( c_xmom_boundary_values   ,gpu_xmom_boundary_values   ,    number_of_xmom_boundary_values  * sizeof(double) );
   copyDeviceToHost( c_ymom_boundary_values   ,gpu_ymom_boundary_values   ,    number_of_ymom_boundary_values  * sizeof(double) );
   */
   copyDeviceToHost( c_stage_explicit_update  ,gpu_stage_explicit_update  ,    number_of_stage_explicit_update* sizeof(double) );
   copyDeviceToHost( c_xmom_explicit_update   ,gpu_xmom_explicit_update   ,    number_of_xmom_explicit_update  * sizeof(double) );
   copyDeviceToHost( c_ymom_explicit_update   ,gpu_ymom_explicit_update   ,    number_of_ymom_explicit_update  * sizeof(double) );
	
   if( iteration == final_iter ) {
   // Free gpu memory
   freeDeviceMemory( gpu_neighbours ); 
   freeDeviceMemory( gpu_neighbour_edges );
   freeDeviceMemory( gpu_normals );
   freeDeviceMemory( gpu_edgelengths);
   freeDeviceMemory( gpu_radii );
   freeDeviceMemory( gpu_areas );
   freeDeviceMemory( gpu_tri_full_flag );
   //freeDeviceMemory( gpu_already_computed_flux);
   freeDeviceMemory( gpu_max_speed_array );
                                                
   freeDeviceMemory( gpu_stage_edge_values );
   freeDeviceMemory( gpu_xmom_edge_values );
   freeDeviceMemory( gpu_ymom_edge_values );
   freeDeviceMemory( gpu_bed_edge_values );
                                                
   freeDeviceMemory( gpu_stage_boundary_values );
   freeDeviceMemory( gpu_xmom_boundary_values );
   freeDeviceMemory( gpu_ymom_boundary_values );
                                                
   freeDeviceMemory( gpu_stage_explicit_update);
   freeDeviceMemory( gpu_xmom_explicit_update );
   freeDeviceMemory( gpu_ymom_explicit_update );
   }

	iteration++;

  Py_DECREF(neighbours);
  Py_DECREF(neighbour_edges);
  Py_DECREF(normals);
  Py_DECREF(edgelengths);
  Py_DECREF(radii);
  Py_DECREF(areas);
  Py_DECREF(tri_full_flag);
  Py_DECREF(already_computed_flux);
  Py_DECREF(max_speed_array);
  Py_DECREF(stage_edge_values);
  Py_DECREF(xmom_edge_values);
  Py_DECREF(ymom_edge_values);
  Py_DECREF(bed_edge_values);
  Py_DECREF(stage_boundary_values);
  Py_DECREF(xmom_boundary_values);
  Py_DECREF(ymom_boundary_values);
  Py_DECREF(stage_explicit_update);
  Py_DECREF(xmom_explicit_update);
  Py_DECREF(ymom_explicit_update);

  
  // Return updated flux timestep
  return Py_BuildValue("d", timestep);
}


//-------------------------------
// Method table for python module
//-------------------------------
static struct PyMethodDef MethodTable[] = {
  /* The cast of the function is necessary since PyCFunction values
   * only take two PyObject* parameters, and rotate() takes
   * three.
   */

  {"compute_fluxes_ext_central_new_gpu", compute_fluxes_ext_central_new_gpu, METH_VARARGS, "Print out"},
  {NULL, NULL}
};

// Module initialisation
void initgpu_python_glue(void){
  Py_InitModule("gpu_python_glue", MethodTable);

  import_array(); // Necessary for handling of NumPY structures
}

