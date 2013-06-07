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
#ifndef OPENACC_H
#define OPENACC_H

#include <stddef.h>

#ifdef OPENACC_API
#  undef OPENACC_API
#endif /* OPENACC_API */

#ifdef _WIN32
#  ifdef OPENACC_BUILD
#    define OPENACC_API __declspec(dllexport)
#  else /* ! OPENACC_BUILD */
//#    pragma comment(lib, "openacc")
#    define OPENACC_API __declspec(dllimport)
#  endif /* OPENACC_BUILD */
#else /* ! _WIN32 */
#  define OPENACC_API
#endif /* _WIN32 */

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * The list of devices types handled by the OPENACC API.
 *
 * @{
 */
typedef enum
{
  acc_device_none     = 0, //< no device
  acc_device_default  = 1, //< default device type
  acc_device_host     = 2, //< host device
  acc_device_not_host = 3, //< not host device
  acc_device_cuda     = 4, //< CUDA device
  acc_device_opencl   = 5  //< OpenCL device
} acc_device_t;
/**
 * @}
 */


/**
 *  Get the number of accelerator devices of the given type.
 *
 * \return the number of accelerator devices of the given type attached to the host.
 *
 * \param devicetype is the kind of device to count.
 */
OPENACC_API
int acc_get_num_devices(acc_device_t devicetype);

/**
 *  Tells the runtime which type of device to use.
 *
 * \param devicetype is the kind of device the runtime must use.
 */
OPENACC_API
void acc_set_device_type(acc_device_t devicetype);

/**
 *  Get the type of the device that will be used to run the next accelerator parallel or kernels region.
 *
 * \return the type of the next used device.
 *
 */
OPENACC_API
acc_device_t acc_get_device_type(void);

/**
 *  Tells the runtime which device to use.
 *
 * \param devicenum is the number of the device the runtime must use (0 means reverting using default behavior).
 * \param devicetype is the kind of device the runtime must use.
 */
OPENACC_API
void acc_set_device_num(int devicenum, acc_device_t devicetype);

/**
 * Get the device number of the specified device type that will be used to execute the next accelerator parallel or kernels region.
 *
 * \return an integer corresponding to the device number of the specified device type.
 * \param devicetype is the kind of device.
 */
OPENACC_API
int acc_get_device_num(acc_device_t devicetype);

/**
 *  Tests for completion of all associated asynchronous activities.
 *
 * \param arg is the integer in the async clause.
 */
OPENACC_API
int acc_async_test(int arg);

/**
 *  Tests for completion of all asynchronous activities.
 */
OPENACC_API
int acc_async_test_all(void);

/**
 *  Waits for completion of all associated asynchronous activities.
 *
 * \param arg is the integer in the async clause.
 */
OPENACC_API
void acc_async_wait(int arg);

/**
 *  Waits for completion of asynchronous activities.
 */
OPENACC_API
void acc_async_wait_all(void);

/**
 *  The acc_init routine tells the runtime to initialize the runtime for that device type.
 *
 * \param devicetype is the kind of device to initialize.
 */
OPENACC_API
void acc_init(acc_device_t devicetype);

/**
 *  The acc_shutdown routine tells the runtime to shut down the connection to the given accelerator device, and free up any runtime resources.
 *
 * \param devicetype is the kind of device to shut down.
 */
OPENACC_API
void acc_shutdown(acc_device_t devicetype);

/**
 *  The acc_on_device routine tells the program whether it is executing on a particular device.
 */
OPENACC_API
int acc_on_device(acc_device_t devicetype);

/**
 *  Allocate memory on the device.
 *
 * \return pointer to the allocated memory on the device.
 * \param size is the size of memory to allocate.
 */
OPENACC_API
void * acc_malloc(size_t size);

/**
 *  free memory on the device.
 *
 * \param ptr is the pointer to the memory to free.
 */
OPENACC_API
void acc_free(void * ptr);

#ifdef __cplusplus
}
#endif

#endif /*OPENACC_H*/
