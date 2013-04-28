
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


    int timestepping_method, flow_algorithm, compute_fluxes_method;
    int skip_initial_step;
    int step;
    double yieldstep, finaltime, duration, epsilon;
    double tmp_timestep;


    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Oddddiiiii", 
                &domain,
                &yieldstep, 
                &finaltime, 
                &duration,
                &epsilon,
                &skip_initial_step,

                &compute_fluxes_method,
                &flow_algorithm,
                &timestepping_method,
                &step
                )) 
    {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }

    static struct domain D;
    
    
    D.timestepping_method = timestepping_method;
    D.flow_algorithm = flow_algorithm;
    D.compute_fluxes_method = compute_fluxes_method;
    

    if ( !step )
    {   
        get_python_domain(&D, domain);
        print_domain_struct(&D);
        //fflush(stdout);
    }
    

    //-------------------------------
    // Start evolve procedure
    //-------------------------------
    tmp_timestep = evolve( &D, 
            yieldstep, finaltime, duration, 
            epsilon, skip_initial_step, step);

    printf(" GLUE: flux_timestep    %lf\n", D.flux_timestep);
    printf(" GLUE: tmp_timestep     %lf\n", tmp_timestep);
    printf(" GLUE: fnltime-epsilon  %lf\n", D.finaltime - epsilon);
    printf(" GLUE: yieldtime        %lf\n", D.yieldtime);
    return Py_BuildValue("d", tmp_timestep);
}



// For testing the python-hmpp linking 
PyObject *hmpp_python_test(PyObject *self, PyObject *args)
{
    PyObject *domain;


    int timestepping_method, flow_algorithm, compute_fluxes_method;
    int skip_initial_step;
    int step;
    double yieldstep, finaltime, duration, epsilon;
    double tmp_timestep;


    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Oddddiiiii", 
                &domain,
                &yieldstep, 
                &finaltime, 
                &duration,
                &epsilon,
                &skip_initial_step,

                &compute_fluxes_method,
                &flow_algorithm,
                &timestepping_method,
                &step
                )) 
    {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }

    static struct domain D;
    
    
    D.timestepping_method = timestepping_method;
    D.flow_algorithm = flow_algorithm;
    D.compute_fluxes_method = compute_fluxes_method;
    

    if ( !step )
    {   
        get_python_domain(&D, domain);
        print_domain_struct(&D);
        //fflush(stdout);
    }
    

    //-------------------------------
    // Start evolve procedure
    //-------------------------------

    // For testing single function
    test_single( &D);
    return Py_BuildValue("");
}



PyObject *hmpp_distribute_to_vertices_and_edges(PyObject *self, PyObject *args)
{
    PyObject *domain;


    int timestepping_method, flow_algorithm, compute_fluxes_method;
    int skip_initial_step;
    int step;
    double yieldstep, finaltime, duration, epsilon;


    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Oddddiiiii", 
                &domain,
                &yieldstep, 
                &finaltime, 
                &duration,
                &epsilon,
                &skip_initial_step,

                &compute_fluxes_method,
                &flow_algorithm,
                &timestepping_method,
                &step
                )) 
    {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }

    static struct domain D;
    
    
    D.timestepping_method = timestepping_method;
    D.flow_algorithm = flow_algorithm;
    D.compute_fluxes_method = compute_fluxes_method;
    

    if ( !step )
    {   
        get_python_domain(&D, domain);
        print_domain_struct(&D);
        //fflush(stdout);
    }
    
    // For testing single function
    //distribute_to_vertices_and_edges(&D);

    _extrapolate_second_order_sw(&D);
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
    {"hmpp_distribute_to_vertices_and_edges", hmpp_distribute_to_vertices_and_edges, METH_VARARGS,"Print out"},
    {NULL, NULL}
};



// Module initialisation
void inithmpp_python_glue(void){
    Py_InitModule("hmpp_python_glue", MethodTable);

    import_array(); // Necessary for handling of NumPY structures
}

