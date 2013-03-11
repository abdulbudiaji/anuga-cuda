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
#ifndef OPENACCI_OPENACCI_C_H
#define OPENACCI_OPENACCI_C_H

#include <stddef.h>

#ifdef OPENACCI_API
#  undef OPENACCI_API
#endif /* OPENACCI_API */

#ifdef _WIN32
#  ifdef OPENACCI_BUILD
#    define OPENACCI_API __declspec(dllexport)
#  else /* ! OPENACCI_BUILD */
#    pragma comment(lib, "openacci")
#    define OPENACCI_API __declspec(dllimport)
#  endif /* OPENACCI_BUILD */
#else /* ! _WIN32 */
#  define OPENACCI_API
#endif /* _WIN32 */

#ifdef __cplusplus
extern "C"
{
#endif

OPENACCI_API
void openacci_call(const char * file_name,
                   int          line_number,
                   const char * function_name);

OPENACCI_API
void openacci_fallback(const char * file_name,
                       int          line_number,
                       const char * function_name);

OPENACCI_API
void openacci_enter_region(const char * file_name,
                           int          line_number,
                           int          region_kind,
                           int          num_args,
                           int          async_mode,
                           int          queue_id);

OPENACCI_API
void openacci_leave_region(const char * file_name,
                           int          line_number);

OPENACCI_API
void openacci_push_data(const char * file_name,
                        int          line_number,
                        const char * variable_name,
                        const void * host_address,
                        size_t       start,
                        size_t       length,
                        size_t       element_size,
                        int          transfer_mode);

OPENACCI_API
void openacci_update_datas(const char *  file_name,
                           int           line_number,
                           int           nb_variables,
                           const char ** variable_names,
                           const void ** host_addresses,
                           size_t     *  starts,
                           size_t     *  lengths,
                           size_t     *  elements_sizes,
                           int        *  update_sides,
                           int           async_mode,
                           int           queue_id);

OPENACCI_API
void openacci_wait(const char * file_name,
                   int          line_number,
                   int          async_mode,
                   int          queue_id);

OPENACCI_API
void * openacci_get_device_pointer(const char * file_name,
                                   int          line_number,
                                   const void * host_address);

#ifdef __cplusplus
}
#endif

#endif
