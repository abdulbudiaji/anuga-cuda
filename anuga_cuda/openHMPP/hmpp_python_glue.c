
#include "Python.h" 
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "numpy_shim.h"


#include "util_ext.h"
#include "sw_domain.h"

PyObject *hmpp_evolve(PyObject *self, PyObject *args)
{
    PyObject * domain, * quantities;
    struct domian D;

    double yieldstep, finaltime, duration;

    
    // Convert Python arguments to C
    if ( !PyArg_ParseTuple( args, "Oddd",
                &domain,
                &yieldstep,
                &finaltime,
                &duration))
    {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }

    get_python_domain(&D, domain);
    print_domain_struct(&D);
}



//-------------------------------
// Method table for python module 
//-------------------------------
static struct PyMethodDef MethodTable[] = {
    {"hmpp_evolve", hmpp_evolve, METH_VARARGS, "Print out"},
    {NULL, NULL}
};


// Module initialization
void inithmpp_python_glue(void){
    Py_InitModule("hmpp_python_glue", MethodTable);

    import_array(); // Necessary for handling of NumPY structures
}
