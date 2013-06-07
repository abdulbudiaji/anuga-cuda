
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
#include "sw_domain.h"

#include "cudafun.h"
#include "evolvefun.h"

//=========================================================================
// Python Glue
//=========================================================================

PyObject *evolve_new_gpu(PyObject *self, PyObject *args) 
{
    PyObject 
        *domain,
        *quantities;


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


    quantities = get_python_object(domain, "quantities");
    if (!quantities) {
        PyErr_SetString(PyExc_RuntimeError, 
                "Can not parse quantities");   
        return NULL;
    }
    

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

    distribute_to_vertices_and_edges(D);
}


//-------------------------------
// Method table for python module
//-------------------------------
static struct PyMethodDef MethodTable[] = {
    /* The cast of the function is necessary since PyCFunction values
     * only take two PyObject* parameters, and rotate() takes
     * three.
     */

    {"evolve_new_gpu", evolve_new_gpu, METH_VARARGS, "Print out"},
    {NULL, NULL}
};

// Module initialisation
void initgpu_python_glue(void){
    Py_InitModule("gpu_python_glue", MethodTable);

    import_array(); // Necessary for handling of NumPY structures
}

