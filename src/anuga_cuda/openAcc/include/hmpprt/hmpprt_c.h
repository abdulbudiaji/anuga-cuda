/*
 * Copyright (C) 2008-2013 CAPS entreprise.  All Rights Reserved.
 * 
 * The source code contained or described herein and all documents  related
 * to the source code ("Material") are owned  by  CAPS  entreprise  or  its
 * suppliers or licensors.
 * 
 * Title to the Material remains with CAPS entreprise or its suppliers  and
 * licensors.  The Material contains trade secrets and proprietary and con-
 * fidential information of CAPS entreprise or its suppliers and licensors.
 * 
 * The Material is protected by the French intellectual property code,  in-
 * tellectual property laws and international treaties.  No part of the Ma-
 * terial may be used, copied, reproduced, modified,  published,  uploaded,
 * posted, transmitted, distributed or disclosed in any  way  without  CAPS
 * entreprise's prior express written permission.
 * 
 * No license under any patent, copyright, trade secret or other  intellec-
 * tual property right is granted to or conferred upon you by disclosure or
 * delivery of the Material, either expressly, by implication,  inducement,
 * estoppel or otherwise.
 * 
 * Any license under such intellectual property rights  must  be  expressed
 * and approved by CAPS entreprise in writing.
 */
#ifndef HMPPRT_C_H
#define HMPPRT_C_H

#include <stddef.h>

#ifdef HMPPRT_API
#  undef HMPPRT_API
#endif /* HMPPRT_API */

#ifdef _WIN32
#  ifdef HMPPRT_BUILD
#    define HMPPRT_API __declspec(dllexport)
#  else /* ! HMPPRT_BUILD */
#    pragma comment(lib, "hmpprt")
#    define HMPPRT_API __declspec(dllimport)
#  endif /* HMPPRT_BUILD */
#else /* ! _WIN32 */
#  define HMPPRT_API
#endif /* _WIN32 */

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *  \defgroup grp_errors     Error Codes
 *  @{
 *  Below are all the possible values for the hmpprt_error_t type used to
 *  represent error status
 *  @}
 *  \defgroup grp_memspace   Memory Spaces
 *  @{
 *  Below are all the possible values for the hmpprt_memspace_t type used to
 *  represent memory spaces
 *  @}
 *  \defgroup grp_targets    Targets
 *  @{
 *  Below are all the possible values for the hmpprt_target_t type used to
 *  represent the known HMPP targets
 *  @}
 *  \defgroup grp_handlers   Object Handlers
 *  @{
 *  Each of the compound types listed below represents an object manipulated
 *  by the HMPP runtime
 *
 *  The meaning of the field \c handle is not documented (this is actually
 *  an address to the real object) but it can be assumed that a value of 0
 *  represents an uninitialized object (e.g. NULL in C or C++)
 *
 *  @}
 *  \defgroup grp_device     Devices
 *  @{
 *  This section presents all functions related to device management
 *  @}
 *  \defgroup grp_grouplet   Grouplets
 *  @{
 *  This section presents all functions related to grouplet management
 *  @}
 *  \defgroup grp_codelet    Codelets
 *  @{
 *  This section presents all functions related to Codelet management
 *  @}
 *  \defgroup grp_arglist    Argument Lists
 *  @{
 *  This section presents all functions related to Argument List management
 *  @}
 *  \defgroup grp_queues     Execution Queues
 *  @{
 *  This section presents all functions related to Queue management
 *  @}
 *  \defgroup grp_data       Data Buffers
 *  @{
 *  This section presents all functions related to Data Buffer management
 *  @}
 *  \defgroup grp_cuda       CUDA Specific Functions
 *  @{
 *  @}
 *  \defgroup grp_opencl     OpenCL Specific Functions
 *  @{
 *  @}
 */

/**
 * \mainpage HMPPRT C API
 *
 * \section intro_sec Introduction
 *
 *  The HMPPRT C API is provided via the header \ref hmpprt_c.h.
 *
 *
 *  Quick links:
 *    - \ref grp_errors
 *    - \ref grp_memspace
 *    - \ref grp_targets
 *    - \ref grp_handlers
 *    - \ref grp_device
 *    - \ref grp_grouplet
 *    - \ref grp_codelet
 *    - \ref grp_arglist
 *    - \ref grp_queues
 *    - \ref grp_data
 *    - \ref grp_cuda
 *    - \ref grp_opencl
 */



/**
 * Describe a HMPP Device
 * \ingroup grp_device
 */
typedef struct HmpprtDevice_t    * hmpprt_device_t ;

/**
 * Describe a HMPP Grouplet
 * \ingroup grp_grouplet
 */
typedef struct HmpprtGrouplet_t  * hmpprt_grouplet_t ;

/**
 * Describe a HMPP Codelet
 * \ingroup grp_codelet
 */
typedef struct HmpprtCodelet_t   * hmpprt_codelet_t ;

/**
 * Describe a HMPP Data buffer
 * \ingroup grp_data
 */
typedef struct HmpprtData_t      * hmpprt_data_t ;

/**
 * Describe a HMPP Argument list
 * \ingroup grp_arglist
 */
typedef struct HmpprtArgList_t   * hmpprt_arglist_t ;

/**
 * Describe a HMPP execution Queue
 * \ingroup grp_queue
 */
typedef struct HmpprtQueue_t     * hmpprt_queue_t ;

/**
 * Describe any address on any device
 *
 * \ingroup grp_handlers
 */
typedef struct HmpprtAddress_t   * hmpprt_address_t ;

/**
 * The list of targets for which a grouplet can be generated.
 *
 * \addtogroup grp_targets
 * @{
 */
// QUESTION: is that a mask?
typedef enum
{
  HMPPRT_ANY_TARGET	= 0, //< Any target (that is not a real target)
  HMPPRT_HOST_TARGET	= 1, //< Host target
  HMPPRT_CUDA_TARGET	= 2, //< CUDA target
  HMPPRT_OPENCL_TARGET	= 4  //< OpenCL target
} hmpprt_target_t ;
/**
 * @}
 */


/**
 * \addtogroup grp_errors
 * @{
 */
typedef enum {

#define HMPP_ERROR_DEF(Klass,Kode,Id) HMPPRT_##Kode##_ERROR = Id,
#include <hmpperr/Error.def>
#undef  HMPP_ERROR_DEF

  HMPPRT_UNKNOWN_ERROR    = 1001,   // An unknown error
  HMPPRT_TIMEOUT_ERROR    = -1,     // Not really an error. indicate a timeout during a wait
  HMPPRT_SUCCESS          = 0       // not an error

} hmpprt_error_t ;
/**
 * @}
 */


/**
 * Memory space
 *
 * \addtogroup grp_memspace
 * @{
 */
typedef enum
{
#define D(id, en, sn) HMPPRT_##en = id,
#include <hmpprt/MemorySpace.def>
#undef D
  HMPPRT_MS_MAX
} hmpprt_memspace_t;
/**
 * @}
 */



/*  ==================================================
 *  ====================  bool  ============
 *  ==================================================
 */

// hmpprt_bool_t, HMPPRT_FALSE and HMPPRT_TRUE are just convenient
// to identify locations where boolean values are expected

typedef int hmpprt_bool_t ;

#define HMPPRT_FALSE 0
#define HMPPRT_TRUE  1


/*  ==================================================
 *  ====================  Error Managment  ===========
 *  ==================================================
 */

/*  ==================================================
 *  ====================  hmpprt_grouplet_t  ============
 *  ==================================================
 */

/**
 *  Load a dynamic library generated by HMPP for a grouplet (or a single codelet)
 *
 *  The filename should include the file extension (e.g. ".so" on Linux or ".dll"
 *  on windows)
 *
 * \param grouplet is the HMPP grouplet to load
 * \param filename is the full filename of the library
 *
 *  \ingroup grp_grouplet
 */
HMPPRT_API
hmpprt_error_t   hmpprt_grouplet_load(hmpprt_grouplet_t * grouplet, const char * filename) ;

/**
 * Similar to hmpprt_grouplet_load but for a grouplet of a specific target.
 *
 * \param filename is the name of the source file used to build the grouplet (with or without the extenstion). Ex: sgemm.c is the source file used to build sgemm_cuda.so, filename should be "sgemm" or "sgemm.c"
 *
 * \param grouplet us the HMPP grouplet to load
 * \param filename is the full filename of the library
 * \param target is the given target
 *
 * \ingroup grp_grouplet
 */
HMPPRT_API
hmpprt_error_t   hmpprt_grouplet_load_t(hmpprt_grouplet_t * grouplet,const char * filename, hmpprt_target_t target) ;

/**
 * Delete a grouplet previously allocated by hmpprt_grouplet_load or
 * hmpprt_grouplet_load_t
 *
 * \ingroup grp_grouplet
 */
HMPPRT_API
hmpprt_error_t   hmpprt_grouplet_delete(hmpprt_grouplet_t grp) ;

/**
 * Provide the filename associated to a grouplet
 *
 * \ingroup grp_grouplet
 */
HMPPRT_API
const char *     hmpprt_grouplet_get_filename(hmpprt_grouplet_t grp) ;

/**
 * Provide the target for which a grouplet was compiled
 *
 * \ingroup grp_grouplet
 */
HMPPRT_API
hmpprt_target_t  hmpprt_grouplet_get_target(hmpprt_grouplet_t grp) ;

/**
 * Provide the number of codelets into a grouplet
 *
 * \ingroup grp_grouplet
 */
HMPPRT_API
int               hmpprt_grouplet_get_codelet_count(hmpprt_grouplet_t grp) ;

/**
 * Provide the codelet of rank n in this grouplet (first codelet is at rank 0) or NULL
 * if an illegal rank is specified
 *
 * Remark: The first codelet has rank 0
 *
 * \ingroup grp_grouplet
 */
HMPPRT_API
hmpprt_codelet_t  hmpprt_grouplet_get_codelet(hmpprt_grouplet_t grp,int n) ;

/**
 * Provide the codelet with the given function name in this grouplet or NULL
 * if an illegal name is specified
 *
 * \ingroup grp_grouplet
 */
HMPPRT_API
hmpprt_codelet_t  hmpprt_grouplet_get_codelet_n(hmpprt_grouplet_t grp,const char *name) ;


/*  ==================================================
 *  ==================== Generic functions     =======
 *  ==================================================
 */

/**
 * Provide a readable name for a supported target
 */
HMPPRT_API
const char *	hmpprt_get_target_name (hmpprt_target_t target);

/*  ==================================================
 *  ==================== Host Device Mgt       =======
 *  ==================================================
 */

/**
 * Provide the unique host device
 */
HMPPRT_API
hmpprt_device_t	hmpprt_get_host_device() ;


/*  ==================================================
 *  ==================== Cuda Device Mgt       =======
 *  ==================================================
 */

/**
 * Fill the array \c devices with up to \c nb CUDA devices.
 *
 * \return the number of CUDA devices filled in \c devices
 *
 * \param nb      is the size of the argument \c devices
 * \param devices is the result
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
int		hmpprt_cuda_get_devices(int nb, hmpprt_device_t* devices) ;

/** Provide the number of cuda devices on the system
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
int		hmpprt_cuda_get_device_count(void)  ;

/**
 *  Get the CUDA device at the specified rank (or NULL)
 *
 *  Numbering starts at 0
 *
 * \param n is the rank of the requested device
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
hmpprt_device_t	hmpprt_cuda_get_device(int n) ;

/**
 * Provide the name of the specified device.
 *
 * The resulting string should not be freed and is insured to remain valid
 * until the \c device is deleted.
 *
 *
 * \ingroup grp_device
 */
HMPPRT_API
const char *	hmpprt_cuda_get_device_name(hmpprt_device_t device) ;

/**
 * Provide the name of the device vendor.
 *
 * The resulting string should not be freed and is insured to remain valid
 * until the \c device is deleted.
 *
 * \ingroup grp_device
 */
HMPPRT_API
const char *	hmpprt_cuda_get_device_vendor(hmpprt_device_t device) ;

/**
 *  Provide the version of the CUDA driver
 *
 * \param device should be a CUDA device
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
int             hmpprt_cuda_get_driver_version(hmpprt_device_t device) ;

/**
 * Provide the compute capability of a CUDA device in *major and *minor.
 *
 * For example, a device with compute capabilities 1.3 would
 * set *major=1 and *minor=3
 *
 * \param device should be a CUDA device
 * \param major  is set to the major part of the device compute capability
 * \param minor  is set to the minor part of the device compute capability
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
void            hmpprt_cuda_get_compute_capability(hmpprt_device_t device, int *major, int *minor) ;

/**
 * Provide the amount of memory available on a CUDA device
 *
 * \remark: 0 is returned for other devices
 *
 * \param device should be a CUDA device
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
size_t     	hmpprt_cuda_get_memory_size(hmpprt_device_t device) ;

/** \return true if memory ECC (error correction code) is enabled for this device
 *
 * \param device should be a CUDA device
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
hmpprt_bool_t 	hmpprt_cuda_has_ecc(hmpprt_device_t device) ;

/** \return true if the CUDA context can be reused for interoperability with CUDA-runtime API.
 *
 * \param device should be a CUDA device
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
hmpprt_bool_t   hmpprt_cuda_has_cuda_interop_capabilities(hmpprt_device_t device) ;

/** Lock the CUDA context for use from CUDA runtime API.
 *
 * \param device should be a CUDA device
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
void            hmpprt_cuda_lock_context(hmpprt_device_t device) ;

/** Unlock the CUDA context for use from CUDA runtime API.
 *
 * \param device should be a CUDA device
 *
 * \ingroup grp_cuda
 */
HMPPRT_API
void            hmpprt_cuda_unlock_context(hmpprt_device_t device) ;

/** Enable/disable the malloc/calloc override by specific CUDA/OPENCL instrinsics.
 *
 * \param replacer intrinsic which must be called instead of malloc/calloc
 *                 (0: host, 1: uva, 2:cl_alloc_host_ptr, 3:cl_use_persistent_mem_amd)

 * \param status on/off
 */
HMPPRT_API
void            hmpprt_overload_mallocs(int replacer, int status) ;

/*  ==================================================
 *  ==================== OpenCL Device Mgmt ==========
 *  ==================================================
 */

/**
 * Fill the array \c devices with up to \c nb OpenCL devices.
 *
 * \param n       is the size of the argument \c devices
 * \param devices is the result
 *
 * \return the number of OpenCL devices filled in \c devices
 *
 * \ingroup grp_opencl
 */
HMPPRT_API
int		hmpprt_opencl_get_devices(int nb, hmpprt_device_t *devices) ;

/**
 *  Get the total number of OpenCL devices
 *
 * \ingroup grp_opencl
 */
HMPPRT_API
int		hmpprt_opencl_get_device_count(void)  ;

/**
 *  Get the OpenCL device at the specified rank (or NULL)
 *
 *  Numbering starts at 0
 *
 * \ingroup grp_opencl
 */
HMPPRT_API
hmpprt_device_t 	hmpprt_opencl_get_device(int rank) ;

/**
 * Provide the name of the specified device.
 *
 * The resulting string should not be freed and is insured to remain valid
 * until the \c device is deleted.
 *
 *
 * \ingroup grp_device
 */
HMPPRT_API
const char *	hmpprt_opencl_get_device_name(hmpprt_device_t device) ;

/**
 * Provide the name of the device vendor.
 *
 * The resulting string should not be freed and is insured to remain valid
 * until the \c device is deleted.
 *
 * \ingroup grp_device
 */
HMPPRT_API
const char *	hmpprt_opencl_get_device_vendor(hmpprt_device_t device) ;

/**
 * \param device should be an OpenCL device
 *
 * \ingroup grp_opencl
 */
HMPPRT_API
const char *   	hmpprt_opencl_get_driver_version(hmpprt_device_t device) ;

/**
 * \param device should be an OpenCL device
 *
 * \ingroup grp_opencl
 */
HMPPRT_API
const char *  	hmpprt_opencl_get_device_version(hmpprt_device_t device) ;

/**
 * \param device should be an OpenCL device
 *
 * \ingroup grp_opencl
 */
HMPPRT_API
size_t     	hmpprt_opencl_get_memory_size(hmpprt_device_t device) ;

/**
 * \param device should be an OpenCL device
 *
 * \ingroup grp_opencl
 */
HMPPRT_API
void *          hmpprt_opencl_get_kernel_queue(hmpprt_device_t device);

/*  ==================================================
 *  ==================== Generic Device Managment ====
 *  ==================================================
 */

/**
 *  Fill the array \c devices with up to \c nb devices of the specified \c target
 *
 *  The returned value is the number devices filled in \c devices
 *
 * \ingroup grp_device
 */
HMPPRT_API
int		hmpprt_get_devices(hmpprt_target_t target, int nb, hmpprt_device_t * devices) ;

/**
 *  Get the total number of devices of the specified \c target
 *
 * \ingroup grp_device
 */
HMPPRT_API
int		hmpprt_get_device_count(hmpprt_target_t target) ;

/**
 *  Get the device of the specified \c target at the specified rank (or NULL)
 *
 *  Numbering starts at 0
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_device_t 	hmpprt_get_device(hmpprt_target_t target, int rank) ;

/*  ==================================================
 *  ==================== Device Mgt            =======
 *  ==================================================
 */

/**
 * Acquire the specified device.
 *
 * The following return codes can be expected
 *   -	HMPPRT_SUCCESS        in case of success
 *   -	HMPPRT_LOCKFAIL_ERROR if this device is already in use.
 *   -	HMPPRT_IO_ERROR       if an error occurred while communicating with the resource manager.
 *   -  HMPPRT_UNKNOWN_ERROR  if an unexpected error occured
 *
 * \see hmpprt_device_try_acquire
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_error_t	hmpprt_device_acquire(hmpprt_device_t device) ;

/**
 * Attempt to acquire the specified device.
 *
 * \return 1 in case of success else 0
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_bool_t	hmpprt_device_try_acquire(hmpprt_device_t device) ;

/**
 * Release a previously acquired device
 *
 * \param device is any previously acquired device
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_error_t	hmpprt_device_release(hmpprt_device_t device);

/**
 * Perform a synchronous call of the specified \c codelet with the specified argument list \c args
 * on the specified \c device
 *
 * \param device  is any device
 * \param codelet is any codelet compatible with that device
 * \param args    is the codelet argument list
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_error_t	hmpprt_device_call(hmpprt_device_t device, hmpprt_codelet_t  codelet, hmpprt_arglist_t args);

/**
 * Provide the default memory space for this device (usually the so-called global memory)
 *
 * \param device  is any device
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_memspace_t hmpprt_device_get_default_memspace(hmpprt_device_t device) ;

/**
 *  Return true if the specified target matches the specified target.
 *
 * \param device is any device
 * \param target is the target to test
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_bool_t   hmpprt_device_match_target(hmpprt_device_t device, hmpprt_target_t target) ;


/**
 * Return true if double floating point precision is supported by this device.
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_bool_t	hmpprt_device_double_support(hmpprt_device_t device) ;

/**
 * Return true if direct access to the host memory is supported by this device.
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_bool_t	hmpprt_device_hostmem_support(hmpprt_device_t device) ;

/**
 * \return the (max) clock rate of the device.
 *
 * The value is given in Mhz.
 *
 * \ingroup grp_device
 */
HMPPRT_API
unsigned int	hmpprt_device_get_max_clock_rate(hmpprt_device_t device) ;

/**
 * \return the unique ID of the device.
 *
 * \ingroup grp_device
 */
HMPPRT_API
unsigned long	hmpprt_device_get_id(hmpprt_device_t device);


/**
 * Allocate memory of the given data
 *
 * \param data the data to allocate
 *
 * \ingroup grp_device
 */
HMPPRT_API
void            hmpprt_device_allocate(hmpprt_data_t data) ;

/**
 * Free memory of the given data.
 *
 * \param data is the data to free
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_error_t  hmpprt_device_free(hmpprt_data_t data) ;

/**
 * Upload data from host to device.
 *
 * \param data is the data to transfer
 * \param host_address is the pointer to the host memory to read from
 * \param offset is the offset of the host_address to read from
 * \param size is the amount of bytes to transfer
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_error_t     hmpprt_device_upload(hmpprt_data_t   data,
                                        void          * host_address,
                                        size_t          offset,
                                        size_t          size ) ;


/**
 * Download data from device to host.
 *
 * \param data is the data to transfer
 * \param host_address is the pointer to the host memory to write to
 * \param offset is the offset of the host_address to write to
 * \param size is the amount of bytes to transfer
 *
 * \ingroup grp_device
 */
HMPPRT_API
hmpprt_error_t     hmpprt_device_download(hmpprt_data_t   data,
                                          void          * host_address,
                                          size_t          offset,
                                          size_t          size ) ;


/*  ==================================================
 *  ====================   hmpprt_codelet_t ==========
 *  ==================================================
 */

/**
 * Provide the grouplet that owns a codelet
 *
 * \param codelet  is the considered codelet
 *
 * \ingroup grp_codelet
 */
HMPPRT_API
hmpprt_grouplet_t hmpprt_codelet_get_grouplet(hmpprt_codelet_t codelet) ;

/**
 * Provide the memory space of an argument specified by its rank
 *
 * \param codelet is the considered codelet
 * \param index   gives the rank of the argument
 *
 * \ingroup grp_codelet
 */
HMPPRT_API
hmpprt_memspace_t hmpprt_codelet_get_arg_memspace_i(hmpprt_codelet_t codelet, int index) ;

// TODO: rename hmpprt_codelet_get_memspace_i ?

/**
 * Provide the memory space of an argument specified by its name.
 *
 * \param codelet  is the codelet
 * \param name     is the name of the argument
 *
 * \ingroup grp_codelet
 */
HMPPRT_API
hmpprt_memspace_t hmpprt_codelet_get_arg_memspace_n(hmpprt_codelet_t codelet, const char *name) ;

/**
 * Provide the name of a codelet
 *
 * \param codelet  is the considered codelet
 *
 * \ingroup grp_codelet
 */
HMPPRT_API
const char *      hmpprt_codelet_get_name(hmpprt_codelet_t codelet) ;

/**
 * Call a codelet on a given device.
 *
 * The codelet and its arguments must be compatible with that device
 *
 * \param codelet  is the codelet to be call
 * \param device   is the device on which the call should occur
 * \param args     is the argument list to use
 *
 * \ingroup grp_codelet
 */
HMPPRT_API
hmpprt_error_t    hmpprt_codelet_call(hmpprt_codelet_t codelet, hmpprt_device_t device, hmpprt_arglist_t args) ;

/**
 * Provide the number of arguments expected by a codelet
 *
 * \param codelet  is any codelet
 *
 * \ingroup grp_codelet
 */
HMPPRT_API
int               hmpprt_codelet_get_arg_count(hmpprt_codelet_t codelet) ;


/*  ==================================================
 *  ==================== Data Mgt  =======
 *  ==================================================
 */

/**
 * Create a device data descriptor of the specified size in the default memory space
 * of the device.
 *
 * The data itself is not yet allocated (see \ref hmpprt_data_allocate_content)
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_data_t     hmpprt_data_new(hmpprt_device_t device, size_t size) ;

/**
 * Create a device data descriptor of the specified size (in bytes) in the
 * specified memory space of the device.
 *
 * The data itself is not yet allocated (see \ref hmpprt_data_allocate_content)
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_data_t     hmpprt_data_new_m(hmpprt_device_t device, size_t size, hmpprt_memspace_t memory_space) ;

HMPPRT_API
void              hmpprt_data_delete(hmpprt_data_t  data) ;

/**
 * Provide the device on which a data descriptor operates
 *
 * \return the device associated with this data.
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_device_t   hmpprt_data_get_device(hmpprt_data_t  data) ;

/**
 * Provide the size of this data, in bytes.
 *
 * \ingroup grp_data
 */
HMPPRT_API
size_t             hmpprt_data_get_size(hmpprt_data_t  data) ;

/**
 * Provide the memory space of this data.
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_memspace_t  hmpprt_data_memspace(hmpprt_data_t  data) ;

/**
 * Get the device address of the data.
 * For advancer users only.
 * The interpretation of the data is device and memory space dependant:
 *  - For CUDA, an address in global memory is a CUDA pointer (not usable to refer objects in the host memory)
 *  - For OpenCL, an address in global memory is a \c cl_mem (a memory buffer descriptor)
 *
 * \ingroup grp_data
 */
HMPPRT_API
void *          hmpprt_data_get_device_address(hmpprt_data_t  data) ;

/**
 * Set the device address to be managed by this data descriptor
 *
 * Remark: This function is for advanced users only.
 *         The recommanded way to associate a device address to a
 *         data descriptor is via an allocation function such
 *         as hmpprt_data_allocate_content()
 *
 * \ingroup grp_data
 */
HMPPRT_API
void            hmpprt_data_set_device_address(hmpprt_data_t  data, void *addr) ;

/**
 * Allocate memory in a data descriptor using the size
 * specified at its creation
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_error_t  hmpprt_data_allocate_content(hmpprt_data_t  data);

/**
 * Free the memory previously allocated using a function
 * such as hmpprt_data_allocate_content()
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_error_t  hmpprt_data_free_content(hmpprt_data_t  data);

/**
 *  Upload a full data region from the given host address.
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_error_t  hmpprt_data_upload(hmpprt_data_t  data, const void * src_host);

/**
 *  Upload part of a data region from the given host address.
 *  Note: In the current implementation, the offset is added to the data device and to src_host.
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_error_t  hmpprt_data_upload_part(hmpprt_data_t  data,  const void * src_host, size_t offset, size_t size );

/**
 *  Downloads a full data region to the given host address.
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_error_t  hmpprt_data_download(hmpprt_data_t  data, void * dest_host);

/**
 *  Downloads part of a data region to the given host address.
 *
 * \ingroup grp_data
 */
HMPPRT_API
hmpprt_error_t  hmpprt_data_download_part(hmpprt_data_t  data, void* dest_host, size_t offset, size_t size );

/** TODO: Partial transfers */

/*  ==================================================
 *  ====================  hmpprt_arglist_t     =======
 *  ==================================================
 */


/**
 * Create an empty argument list
 *
 * \ingroup grp_arglist
 */
HMPPRT_API
hmpprt_arglist_t hmpprt_arglist_new(void);

/**
 * Create an exact copy of an argument list
 *
 * \param args  is the cloned list
 * \return      the clone list
 *
 * \ingroup grp_arglist
 */
HMPPRT_API
hmpprt_arglist_t hmpprt_arglist_clone(hmpprt_arglist_t args);

/**
 * Provide the number of elements currently stored in the argument list
 *
 * \param args  is an argument list

 * \ingroup grp_arglist
 */
HMPPRT_API
int             hmpprt_arglist_size(hmpprt_arglist_t args);

/**
 * Append the given hmpprt_data_t object to the argument list.
 *
 * \param args  is the argument list to be extended
 * \param data  is the Data buffer

 * \ingroup grp_arglist
 */
HMPPRT_API
void            hmpprt_arglist_add_data(hmpprt_arglist_t args, hmpprt_data_t  data);

/**
 * Append a copy of the \c size bytes starting at data to the argument list.
 *
 * \param args    is the argument list to be extended
 * \param data    is any object whose data must be passed
 * \param size    is the size (in bytes)
 *
 * \ingroup grp_arglist
 */
HMPPRT_API
void            hmpprt_arglist_add_copy(hmpprt_arglist_t args, void * data, size_t size);

/**
 * Append a \c int value to the argument list.
 *
 * This is equivalent to hmpprt_arglist_add_copy(args,&value,sizeof(int))
 *
 * \param args   is the argument list to be extended
 * \param value  is the value to be added
 *
 * \ingroup grp_arglist
 */
HMPPRT_API
void            hmpprt_arglist_add_int(hmpprt_arglist_t args, int value);

/**
 * Append a \c long value to the argument list.
 *
 * This is equivalent to hmpprt_arglist_add_copy(args,&value,sizeof(long))
 *
 * \param args   is the argument list to be extended
 * \param value  is the value to be added
 *
 * \ingroup grp_arglist
 */
HMPPRT_API
void            hmpprt_arglist_add_long(hmpprt_arglist_t args, long value);

/**
 * Append a \c float value to the argument list.
 *
 * This is equivalent to hmpprt_arglist_add_copy(args,&value,sizeof(float))
 *
 * \param args   is the argument list to be extended
 * \param value  is the value to be added
 *
 * \ingroup grp_arglist
 */
HMPPRT_API
void            hmpprt_arglist_add_float(hmpprt_arglist_t args, float value);

/**
 * Append a \c double value to the argument list.
 *
 * This is equivalent to hmpprt_arglist_add_copy(args,&value,sizeof(double))
 *
 * \param args   is the argument list to be extended
 * \param value  is the value to be added
 *
 * \ingroup grp_arglist
 */
HMPPRT_API
void            hmpprt_arglist_add_double(hmpprt_arglist_t args, double value);


/**
 * A special value used as marker to help debug faulty formats
 * (see hmpprt_arglist_fadd)
 */
#define HMPPRT_ARG_CHECK 12345678

/**
 *  Append multiple arguments each specified by a single character of the \c format.
 *
 *  The supported characters are
 *    - '@' to append an argument of type <tt>hmpprt_data_t</tt>
 *    - 'i' to append an argument of type <tt>int</tt> or equivalent
 *    - 'l' to append an argument of type <tt>long</tt> or equivalent
 *    - 's' to append an argument of type <tt>size_t</tt> or equivalent
 *    - 'f' to append an argument of type <tt>float</tt> (see remark below)
 *    - 'd' to append an argument of type <tt>double</tt>
 *    - 'm' to append a copy of a memory area specified by its start address of type <tt>void*</tt>
 *      and its size in byte of type <tt>size_t</tt>
 *    - '#' to insert a marker of value HMPPRT_ARG_CHECK (for debugging purpose)
 *
 *  Example: The following sequence of calls
 *
 * \verbatim
 *     int             n  ;
 *     size_t          len ;
 *     hmpprt_arglist_t  args ;
 *     hmpprt_data_t  T_data  ;
 *     ...
 *     hmpprt_arglist_add_int(args,n)  ;
 *     hmpprt_arglist_add_data(args,T_data) ;
 *     hmpprt_arglist_add_copy(args,len,sizeof(len)) ;
 * \endverbatim
 *
 *  is equivalent to
 *
 *
 * \verbatim
 *     int             n  ;
 *     size_t          len ;
 *     hmpprt_arglist_t  args ;
 *     hmpprt_data_t  T_data  ;
 *     ...
 *     hmpprt_arglist_add_multi(args,"i@m",
 *                              n,
 *                              T_data,
 *                              &len,sizeof(size_t)
 *                             )  ;
 * \endverbatim
 *
 * \ingroup grp_arglist
 */
HMPPRT_API
void            hmpprt_arglist_add_multi(hmpprt_arglist_t args, const char *format, ... );

#if 0

// Those functions are not yet implemented

void            hmpprt_arglist_set_data(hmpprt_arglist_t args, int rank, hmpprt_data_t  data);
void            hmpprt_arglist_set_value(hmpprt_arglist_t args, int rank, void * data, size_t size);
void            hmpprt_arglist_set_int(hmpprt_arglist_t args, int rank, int value);
void            hmpprt_arglist_set_long(hmpprt_arglist_t args, int rank, long value);
void            hmpprt_arglist_set_float(hmpprt_arglist_t args, int rank, float value);
void            hmpprt_arglist_set_double(hmpprt_arglist_t args, int rank, double value);

#endif

/*  ==================================================
 *  ==================== HmpprtQueue  =======
 *  ==================================================
 */

/**
 * Create a new execution queue
 *
 * \ingroup grp_queues
 */
HMPPRT_API
hmpprt_queue_t    hmpprt_queue_new(void) ;

/**
 * Delete an execution queue previously created with hmpprt_queue_new
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_queue_delete(hmpprt_queue_t queue) ;

/**
 * Launch the execution of queue operations and returns immediately.
 * Once the queue is running, it is still possible to enqueue an operation.
 *
 * \param queue   is any queue
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_queue_start(hmpprt_queue_t queue) ;

/**
 * Wait until every operation in the queue is finished
 *
 * \param queue   is any queue
 *
 * \ingroup grp_queues
 */
HMPPRT_API
hmpprt_error_t  hmpprt_queue_wait(hmpprt_queue_t queue) ;

/**
 * Wait until every operation in the queue is finished or until
 * the specified timeout
 *
 * \param queue      is any queue
 * \param timeout_ms is in miliseconds (use 0 to wait for ever)
 *
 * \ingroup grp_queues
 */
HMPPRT_API
hmpprt_error_t  hmpprt_queue_wait_t(hmpprt_queue_t queue, size_t timeout_ms) ;

/**
 * Synchronous execution of the queue. This is an efficient shortcut to
 *
 *     hmpprt_queue_start(queue) ;
 *     hmpprt_queue_wait(queue,0) ;
 *
 * \param queue   is any queue
 *
 * \ingroup grp_queues
 */
HMPPRT_API
hmpprt_error_t  hmpprt_queue_execute(hmpprt_queue_t queue);

/**
 * Clear the content of this queue, removing every operations without executing them.
 *
 * \param queue   is any queue
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_queue_clear(hmpprt_queue_t queue) ;



/**
 * Enqueue the execution of hmpprt_device_acquire(dev)
 *
 * \param queue   is any queue
 * \param device  is any device
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_acquire	(hmpprt_queue_t queue, hmpprt_device_t  device);

/**
 * Enqueue the execution of hmpprt_device_release(dev)
 *
 * \param queue   is any queue
 * \param device  is any device
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_release	(hmpprt_queue_t queue, hmpprt_device_t  device);

/**
 * Enqueue the execution of hmpprt_device_call(dev,cdt,args)
 *
 * \param queue   is any queue
 * \param dev     is any device
 * \param cdt     is any codelet
 * \param args    is an argument list
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_call	(hmpprt_queue_t queue, hmpprt_device_t  dev, hmpprt_codelet_t cdt, hmpprt_arglist_t args);

/**
 *  Enqueue the execution of hmpprt_data_allocate(data)
 *
 * \param queue   is any queue
 * \param data    is any data buffer descriptor
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_allocate	(hmpprt_queue_t queue, hmpprt_data_t  data);

/**
 *  Enqueue the execution of hmpprt_data_free(data)
 *
 * \param queue   is any queue
 * \param data    is any data buffer descriptor
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_free     (hmpprt_queue_t queue, hmpprt_data_t  data);

/**
 *  Enqueue the execution of hmpprt_data_upload(data,host_address)
 *
 * \param queue        is any queue
 * \param data         is any data buffer descriptor
 * \param host_address is the host address to read from
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_upload   (hmpprt_queue_t queue, hmpprt_data_t  data, const void * host_address);

/**
 *  Enqueue the execution of hmpprt_data_upload_part(data,host_address,offset,size)
 *
 * \param queue        is any queue
 * \param data         is any data buffer descriptor
 * \param host_address is any host_address
 * \param offset       is any offset
 * \param size         is any size
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_upload_part   (hmpprt_queue_t queue, hmpprt_data_t  data, const void * host_address, size_t offset, size_t size);

/**
 *  Enqueue the execution of hmpprt_data_download(data,host_address)
 *
 * \param queue        is any queue
 * \param data         is any data buffer descriptor
 * \param host_address is any host_address
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_download (hmpprt_queue_t queue, hmpprt_data_t  data, void * host_address);

/**
 *  Enqueue the execution of hmpprt_data_download(data,host_address,offset,size)
 *
 * \param queue        is any queue
 * \param data         is any data buffer descriptor
 * \param host_address is any host_address
 * \param offset       is any offset
 * \param size         is any size
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_download_part (hmpprt_queue_t queue, hmpprt_data_t  data, void * host_address, size_t offset, size_t size);

/**
 *  Enqueue the execution of hmpprt_queue_execute(sub_queue)
 *
 * \param queue     is any queue
 * \param sub_queue is any sub-queue
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_execute  (hmpprt_queue_t queue, hmpprt_queue_t  sub_queue);

/**
 *  Enqueue the execution of hmpprt_queue_start(sub_queue)
 *
 * \param queue   is any queue
 * \param sub_queue is any sub-queue
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_start    (hmpprt_queue_t queue, hmpprt_queue_t  sub_queue);

/**
 *  Enqueue the execution of hmpprt_queue_start(sub_queue,timeout_ms)
 *
 * \param queue   is any queue
 * \param sub_queue is any sub-queue
 * \param timeout_ms is the time-out (in ms)
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_wait     (hmpprt_queue_t queue, hmpprt_queue_t  sub_queue, size_t timeout_ms);

/**
 *  Enqueue the execution of hmpprt_queue_clear(sub_queue)
 *
 * \param queue     is any queue
 * \param sub_queue is sub-queue
 *
 * \ingroup grp_queues
 */
HMPPRT_API
void            hmpprt_enqueue_clear    (hmpprt_queue_t queue, hmpprt_queue_t  sub_queue);


/**
 *  Return the current time in microseconds
 */
HMPPRT_API
unsigned long hmpprt_get_time_microseconds();

/**
 *  Return the current time in cycles of current CPU
 */
HMPPRT_API
unsigned long hmpprt_get_time_cycles();

#ifdef __cplusplus
} // of extern C
#endif

#endif
