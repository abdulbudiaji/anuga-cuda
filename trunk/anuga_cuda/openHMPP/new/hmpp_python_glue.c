
// Python - C extension module for hmpp_domain.py
//


#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "numpy_shim.h"
#include "util_ext.h"
//#include "sw_domain.h"

#include "hmpp_fun.h"
#include "sw_domain_fun.h"



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

    evolve(D, yieldstep, finaltime, duration, skip_initial_step);
    return Py_BuildValue("");
}



// For testing the python-hmpp linking 
PyObject *hmpp_python_test()
{
    test_call();
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
    {"hmpp_python_test", hmpp_python_test, METH_VARARGS, "Print out"},
    {NULL, NULL}
};



// Module initialisation
void inithmpp_python_glue(void){
    Py_InitModule("hmpp_python_glue", MethodTable);

    import_array(); // Necessary for handling of NumPY structures
}

