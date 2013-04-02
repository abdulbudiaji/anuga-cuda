
// Python - C extension module for hmpp_domain.py
//


#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "numpy_shim.h"
#include "util_ext.h"
#include "sw_domain.h"

#include "gravity.c"


//=========================================================================
// Python Glue
//=========================================================================

PyObject *hmpp_evolve(PyObject *self, PyObject *args) 
{
    PyObject *domain;


    double yieldstep, finaltime, duration;
    int skip_initial_step;


    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Odddb", 
                &domain,
                &yieldstep, 
                &finaltime, 
                &duration,
                &skip_initial_step 
                )) 
    {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }

    struct domain D;
    get_python_domain(&D, domain);
    print_domain_struct(&D);
    //fflush(stdout);

    
    

    //-------------------------------
    // Start evolve procedure
    //-------------------------------

    printf("beta_w %lf\n", D.beta_w);
    if ( D.beta_w < 0 || D.beta_w > 2.0 )
    {
        PyErr_SetString(PyExc_RuntimeError,
                "Attribute self.beta_w must be in the interval [0, 2]\n");
        return NULL;
    }
    
    gravity_call(
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

    return Py_BuildValue("");

}


//-------------------------------
// Method table for python module
//-------------------------------
static struct PyMethodDef MethodTable[] = {
    /* The cast of the function is necessary since PyCFunction values
     * only take two PyObject* parameters, and rotate() takes
     * three.
     */

    {"hmpp_evolve", hmpp_evolve, METH_VARARGS, "Print out"},
    {NULL, NULL}
};

// Module initialisation
void inithmpp_python_glue(void){
    Py_InitModule("hmpp_python_glue", MethodTable);

    import_array(); // Necessary for handling of NumPY structures
}

