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
#ifndef OPENACCI_OPENACCI_FORTRAN_H
#define OPENACCI_OPENACCI_FORTRAN_H

#include <stddef.h>

#ifdef OPENACCI_API
#  undef OPENACCI_API
#endif /* OPENACCI_API */

#ifdef _WIN32
#  ifdef OPENACCI_BUILD
#    define OPENACCI_API __declspec(dllexport)
#  else /* ! OPENACCI_BUILD */
#    define OPENACCI_API __declspec(dllimport)
#  endif /* OPENACCI_BUILD */
#else /* ! _WIN32 */
#  define OPENACCI_API
#endif /* _WIN32 */

#ifdef _WIN32
#  ifdef ABSOFT
#    ifdef _WIN64
#    else // ! _WIN64
#    endif
#  else // ! ABSOFT
#define openacci_call_ OPENACCI_CALL
#define openacci_fallback_ OPENACCI_FALLBACK
#define openacci_enter_global_region_ OPENACCI_ENTER_GLOBAL_REGION
#define openacci_enter_region_ OPENACCI_ENTER_REGION
#define openacci_leave_region_ OPENACCI_LEAVE_REGION
#define openacci_push_global_data_ OPENACCI_PUSH_GLOBAL_DATA
#define openacci_is_initialized_ OPENACCI_IS_INITIALIZED
#define openacci_record_module_ OPENACCI_RECORD_MODULE
#define openacci_push_data_ OPENACCI_PUSH_DATA
#define openacci_update_datas_ OPENACCI_UPDATE_DATAS
#define openacci_allocate_ OPENACCI_ALLOCATE
#define openacci_wait_ OPENACCI_WAIT
#define openacci_get_device_pointer_ OPENACCI_GET_DEVICE_POINTER
#  endif // ! ABSOFT
#endif

#ifdef __cplusplus
extern "C"
{
#endif

OPENACCI_API
void openacci_call_(const char * file_name,
                    const int  & line_number,
                    const char * function_name,
                    int          len_file_name,
                    int          len_function_name);

OPENACCI_API
void openacci_fallback_(const char * file_name,
                        const int  & line_number,
                        const char * function_name,
                        int          len_file_name,
                        int          len_function_name);

OPENACCI_API
void openacci_enter_global_region_(const char * file_name,
                                   const int  & line_number,
                                   const int  & region_kind,
                                   const int  & num_args,
                                   const int  & async_mode,
                                   const int  & queue_id,
                                   int          len_file_name);

OPENACCI_API
void openacci_enter_region_(const char * file_name,
                            const int  & line_number,
                            const int  & region_kind,
                            const int  & num_args,
                            const int  & async_mode,
                            const int  & queue_id,
                            int          len_file_name);

OPENACCI_API
void openacci_leave_region_(const char * file_name,
                            const int  & line_number,
                            int          len_file_name);

OPENACCI_API
bool openacci_is_initialized_(const void   * module_address);

OPENACCI_API
void openacci_record_module_(const void   * module_address);

OPENACCI_API
void openacci_push_global_data_(const char   * file_name,
                                const int    & line_number,
                                const char   * variable_name,
                                const void   * host_address,
                                const size_t & start,
                                const size_t & length,
                                const size_t & element_size,
                                const int    & transfer_mode,
                                int            len_file_name,
                                int            len_variable_name);

OPENACCI_API
void openacci_push_data_(const char   * file_name,
                         const int    & line_number,
                         const char   * variable_name,
                         const void   * host_address,
                         const size_t & start,
                         const size_t & length,
                         const size_t & element_size,
                         const int    & transfer_mode,
                         int            len_file_name,
                         int            len_variable_name);

OPENACCI_API
void openacci_update_datas_(const char   * file_name,
                            const int    & line_number,
                            const int    & nb_variables,
                            const char   * variable_names,
                            const size_t * starts,
                            const size_t * lengths,
                            const size_t * elements_sizes,
                            const int    * update_sides,
                            const int    & async_mode,
                            const int    * queue_id,
                            .../*[const void   * host_addresses,]*
                            int            len_file_name,
                            int            len_variable_names*/);

OPENACCI_API
void openacci_allocate_(const char * file_name,
                        const int  & line_number,
                        const void * host_address,
                        const int  * num_argument,
                        ...);

OPENACCI_API
void openacci_wait_(const char * file_name,
                    const int  & line_number,
                    const int  & async_mode,
                    const int  & queue_id,
                    int          len_file_name);

OPENACCI_API
void openacci_get_device_pointer_(const char * file_name,
                                  int        & line_number,
                                  void      ** device_address,
                                  void       * host_address,
                                  int          len_file_name);

#ifdef __cplusplus
}
#endif

#endif
