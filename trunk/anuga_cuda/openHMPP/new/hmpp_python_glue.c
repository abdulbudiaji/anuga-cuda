
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

#ifdef DEBUG
#define DEBUG_LOG_PAR(a, b) printf(a, b)
#else
#define DEBUG_LOG_PAR(a, b) 
#endif



//=========================================================================
// Python Glue
//=========================================================================

PyObject *hmpp_evolve(PyObject *self, PyObject *args) 
{
    PyObject *domain;


    int timestepping_method, flow_algorithm, compute_fluxes_method;
    int skip_initial_step;
    int step;
    long Nid;
    double yieldstep, finaltime, duration, epsilon;
    double tmp_timestep;


    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Oddddiiiili", 
                &domain,
                &yieldstep, 
                &finaltime, 
                &duration,
                &epsilon,
                &skip_initial_step,

                &compute_fluxes_method,
                &flow_algorithm,
                &timestepping_method,
                &Nid,
                &step
                )) 
    {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }

    static struct domain D;
    
    D.boundary_number = Nid;
    
    D.timestepping_method = timestepping_method;
    D.flow_algorithm = flow_algorithm;
    D.compute_fluxes_method = compute_fluxes_method;
    

    if ( !step )
    {   
        get_python_domain(&D, domain);
        //print_domain_struct(&D);
        //fflush(stdout);
    }
    

    //-------------------------------
    // Start evolve procedure
    //-------------------------------
    tmp_timestep = evolve( &D, 
            yieldstep, finaltime, duration, 
            epsilon, skip_initial_step, step);
/*
    update_boundary(&D);
*/
    DEBUG_LOG_PAR(" GLUE: flux_timestep    %lf\n", D.flux_timestep);
    DEBUG_LOG_PAR(" GLUE: tmp_timestep     %lf\n", tmp_timestep);
    DEBUG_LOG_PAR(" GLUE: fnltime-epsilon  %lf\n", D.finaltime - epsilon);
    DEBUG_LOG_PAR(" GLUE: yieldtime        %lf\n", D.yieldtime);
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
    _distribute_to_vertices_and_edges(&D);
    //update_ghosts(&D);
    //_distribute_to_vertices_and_edges(&D);
    //update_boundary(&D);
    //update_extrema(&D);

    //_extrapolate_second_order_sw(&D);
    return Py_BuildValue("");
}



PyObject *hmpp_extrapolate_second_order_sw(PyObject *self, PyObject *args)
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
        //print_domain_struct(&D);
        //fflush(stdout);
    }
    
    // For testing single function
    _extrapolate_second_order_sw(&D);
    return Py_BuildValue("");
}



PyObject *hmpp_extrapolate_second_order_and_limit_by_vertex(
        PyObject *self, PyObject *args)
{
    PyObject *domain;


    int timestepping_method, flow_algorithm, compute_fluxes_method;
    int skip_initial_step;
    int step;
    long Nid;
    double yieldstep, finaltime, duration, epsilon;


    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Oddddiiiili", 
                &domain,
                &yieldstep, 
                &finaltime, 
                &duration,
                &epsilon,
                &skip_initial_step,

                &compute_fluxes_method,
                &flow_algorithm,
                &timestepping_method,
                &Nid,
                &step
                )) 
    {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }

    static struct domain D;
    
    D.boundary_number = Nid;
    
    D.timestepping_method = timestepping_method;
    D.flow_algorithm = flow_algorithm;
    D.compute_fluxes_method = compute_fluxes_method;
    

    if ( !step )
    {   
        get_python_domain(&D, domain);
        //fflush(stdout);
    }
    
    // For testing single function
    test_extrapolate_second_order_and_limit_by_vertex(&D);
    return Py_BuildValue("");
}

PyObject *hmpp_extrapolate_second_order_and_limit_by_vertex_normal(
        PyObject *self, PyObject *args)
{
    PyObject *domain;


    int timestepping_method, flow_algorithm, compute_fluxes_method;
    int skip_initial_step;
    int step;
    long Nid;
    double yieldstep, finaltime, duration, epsilon;


    // Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Oddddiiiili", 
                &domain,
                &yieldstep, 
                &finaltime, 
                &duration,
                &epsilon,
                &skip_initial_step,

                &compute_fluxes_method,
                &flow_algorithm,
                &timestepping_method,
                &Nid,
                &step
                )) 
    {
        report_python_error(AT, "could not parse input arguments");
        return NULL;
    }
    static struct domain D;
    
    D.boundary_number = Nid;
    
    D.timestepping_method = timestepping_method;
    D.flow_algorithm = flow_algorithm;
    D.compute_fluxes_method = compute_fluxes_method;
    

    if ( !step )
    {   
        get_python_domain(&D, domain);
        //fflush(stdout);
    }
    
    // For testing single function
    test_extrapolate_second_order_and_limit_by_vertex_normal(&D);
    return Py_BuildValue("");
}



PyObject *hmpp_compute_fluxes(PyObject *self, PyObject *args)
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
        //print_domain_struct(&D);
        //fflush(stdout);
    }
    
    // For testing single function
    compute_fluxes(&D);
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
    {"hmpp_extrapolate_second_order_sw", hmpp_extrapolate_second_order_sw, METH_VARARGS, "Print out"},
    {"hmpp_compute_fluxes", hmpp_compute_fluxes, METH_VARARGS, "Print out"},
    {"hmpp_extrapolate_second_order_and_limit_by_vertex", hmpp_extrapolate_second_order_and_limit_by_vertex, METH_VARARGS, "Print out"},
    {"hmpp_extrapolate_second_order_and_limit_by_vertex_normal", hmpp_extrapolate_second_order_and_limit_by_vertex_normal, METH_VARARGS, "Print out"},
    {NULL, NULL}
};



// Module initialisation
void inithmpp_python_glue(void){
    Py_InitModule("hmpp_python_glue", MethodTable);

    import_array(); // Necessary for handling of NumPY structures
}

