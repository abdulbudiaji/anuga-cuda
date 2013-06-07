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
#ifndef HMPPRTI_HMPPRTI_FORTRAN_H
#define HMPPRTI_HMPPRTI_FORTRAN_H

#include <stddef.h>

#ifdef HMPPRTI_API
#  undef HMPPRTI_API
#endif /* HMPPRTI_API */

#ifdef _WIN32
#  ifdef HMPPRTI_BUILD
#    define HMPPRTI_API __declspec(dllexport)
#  else /* ! HMPPRTI_BUILD */
#    define HMPPRTI_API __declspec(dllimport)
#  endif /* HMPPRTI_BUILD */
#else /* ! _WIN32 */
#  define HMPPRTI_API
#endif /* _WIN32 */


#ifdef _WIN32
#  ifdef ABSOFT
#    ifdef _WIN64
#    else // ! _WIN64
#    endif
#  else // ! ABSOFT
#define hmpprti_get_error_ HMPPRTI_GET_ERROR
#define hmpprti_set_error_ HMPPRTI_SET_ERROR
#define hmpprti_lookup_grouplet_ HMPPRTI_LOOKUP_GROUPLET
#define hmpprti_set_property_ HMPPRTI_SET_PROPERTY
#define hmpprti_is_initialized_ HMPPRTI_IS_INITIALIZED
#define hmpprti_set_initialized_ HMPPRTI_SET_INITIALIZED
#define hmpprti_append_target_ HMPPRTI_APPEND_TARGET
#define hmpprti_append_buffer_ HMPPRTI_APPEND_BUFFER
#define hmpprti_append_codelet_ HMPPRTI_APPEND_CODELET
#define hmpprti_append_buffer_to_codelet_ HMPPRTI_APPEND_BUFFER_TO_CODELET
#define hmpprti_append_buffer_to_codelet_hmpp2_ HMPPRTI_APPEND_BUFFER_TO_CODELET_HMPP2
#define hmpprti_push_acquire_ HMPPRTI_PUSH_ACQUIRE
#define hmpprti_push_release_ HMPPRTI_PUSH_RELEASE
#define hmpprti_push_free_all_ HMPPRTI_PUSH_FREE_ALL
#define hmpprti_push_allocate_ HMPPRTI_PUSH_ALLOCATE
#define hmpprti_push_free_ HMPPRTI_PUSH_FREE
#define hmpprti_push_bind_ HMPPRTI_PUSH_BIND
#define hmpprti_push_load_ HMPPRTI_PUSH_LOAD
#define hmpprti_push_wait_load_ HMPPRTI_PUSH_WAIT_LOAD
#define hmpprti_push_store_ HMPPRTI_PUSH_STORE
#define hmpprti_push_wait_store_ HMPPRTI_PUSH_WAIT_STORE
#define hmpprti_push_callsite_ HMPPRTI_PUSH_CALLSITE
#define hmpprti_push_callsite_hmpp2_ HMPPRTI_PUSH_CALLSITE_HMPP2
#define hmpprti_push_standalone_callsite_ HMPPRTI_PUSH_STANDALONE_CALLSITE
#define hmpprti_push_synchronize_ HMPPRTI_PUSH_SYNCHRONIZE
#define hmpprti_push_host_read_ HMPPRTI_PUSH_HOST_READ
#define hmpprti_push_host_write_ HMPPRTI_PUSH_HOST_WRITE
#define hmpprti_push_allocate_mirror_ HMPPRTI_PUSH_ALLOCATE_MIRROR
#define hmpprti_push_free_mirror_ HMPPRTI_PUSH_FREE_MIRROR
#define hmpprti_push_load_mirror_ HMPPRTI_PUSH_LOAD_MIRROR
#define hmpprti_push_wait_load_mirror_ HMPPRTI_PUSH_WAIT_LOAD_MIRROR
#define hmpprti_push_store_mirror_ HMPPRTI_PUSH_STORE_MIRROR
#define hmpprti_push_wait_store_mirror_ HMPPRTI_PUSH_WAIT_STORE_MIRROR
#define hmpprti_run_directive_ HMPPRTI_RUN_DIRECTIVE
#define hmpprti_enqueue_directive_ HMPPRTI_ENQUEUE_DIRECTIVE
#define hmpprti_run_queue_ HMPPRTI_RUN_QUEUE
#  endif // ! ABSOFT
#endif // WIN32

#ifdef __cplusplus
extern "C"
{
#endif

HMPPRTI_API
int hmpprti_get_error_(const char * name,
                       int          len_name);

HMPPRTI_API
void hmpprti_set_error_(const char * name,
                        const int  & error,
                        int          len_name);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_lookup_grouplet_(const char * grouplet_name,
                                                   int          len_grouplet_name);

HMPPRTI_API
void hmpprti_set_property_(const char * name,
                           const char * value,
                           int          len_name,
                           int          len_value);

HMPPRTI_API
int hmpprti_is_initialized_(const struct hmpprti_grouplet * const & grouplet);

HMPPRTI_API
void hmpprti_set_initialized_(struct hmpprti_grouplet * const & grouplet,
                              const int                       & initialized);

HMPPRTI_API
void hmpprti_append_target_(struct hmpprti_grouplet * const & grouplet,
                            const char                      * target_name,
                            int                               len_target_name);

HMPPRTI_API
void hmpprti_append_buffer_(struct hmpprti_grouplet * const & grouplet,
                            const char                      * buffer_name,
                            const size_t                    & element_size,
                            const int                       & dimension_count,
                            const int                       & dynamic,
                            const int                       & lazy_transfer_at_callsite,
                            const int                       & load_once_between_two_allocates,
                            int                               len_buffer_name);

HMPPRTI_API
void hmpprti_append_codelet_(struct hmpprti_grouplet * const & grouplet,
                             const char                      * codelet_name,
                             const char                      * function_name,
                             int                               len_codelet_name,
                             int                               len_function_name);

HMPPRTI_API
void hmpprti_append_buffer_to_codelet_(struct hmpprti_grouplet * const & grouplet,
                                       const int                       & codelet_index,
                                       const int                       & buffer_index,
                                       const int                       & load_before_callsite,
                                       const int                       & store_after_callsite);

HMPPRTI_API
void hmpprti_append_buffer_to_codelet_hmpp2_(struct hmpprti_grouplet * const & grouplet,
                                             const int                       & codelet_index,
                                             const int                       & buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_acquire_(struct hmpprti_grouplet * const & grouplet,
                                                const int                       & default_device,
                                                const int                       & device,
                                                const int                       & exclusive);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_release_(struct hmpprti_grouplet * const & grouplet,
                                                const int                       & default_device,
                                                const int                       & device);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_free_all_(struct hmpprti_grouplet * const & grouplet);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_allocate_(struct hmpprti_grouplet * const & grouplet,
                                                 const int                       & buffer_index,
                                                 const int                       & ignore_unbound,
                                                 const int                       & default_device,
                                                 const int                       & device,
                                                 const int                       * dimension_count,
                                                 ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_free_(struct hmpprti_grouplet * const & grouplet,
                                             const int                       & buffer_index,
                                             const int                       & ignore_unbound);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_bind_(struct hmpprti_grouplet * const & grouplet,
                                             const int                       & buffer_index,
                                             const char                      * address_string,
                                             void                            * host_address,
                                             int                               len_address_string);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_load_(struct hmpprti_grouplet * const & grouplet,
                                             const int                       & buffer_index,
                                             const int                       & async,
                                             const int                       * slice_count,
                                             ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_wait_load_(struct hmpprti_grouplet * const & grouplet,
                                                  const int                       & buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_store_(struct hmpprti_grouplet * const & grouplet,
                                              const int                       & buffer_index,
                                              const int                       & async,
                                              const int                       * slice_count,
                                              ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_wait_store_(struct hmpprti_grouplet * const & grouplet,
                                                   const int                       & buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_callsite_(struct hmpprti_grouplet * const & grouplet,
                                                 const int                       & codelet_index,
                                                 const int                       & async);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_callsite_hmpp2_(struct hmpprti_grouplet * const & grouplet,
                                                       const int                       & codelet_index,
                                                       const int                       & async,
                                                       const int                       * argument_count,
                                                       ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_standalone_callsite_(struct hmpprti_grouplet * const & grouplet,
                                                            const char                      * function_name,
                                                            const int                       & argument_count,
                                                            int                               len_function_name);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_synchronize_(struct hmpprti_grouplet * const & grouplet,
                                                    const int                       & codelet_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_host_read_(struct hmpprti_grouplet * const & grouplet,
                                                  const int                       & buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_host_write_(struct hmpprti_grouplet * const & grouplet,
                                                   const int                       & buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_allocate_mirror_(struct hmpprti_grouplet * const & grouplet,
                                                        const char                      * address_string,
                                                        const int                       & default_device,
                                                        const int                       & device,
                                                        const size_t                    & element_size,
                                                        const int                       * dimension_count,
                                                        ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_free_mirror_(struct hmpprti_grouplet * const & grouplet,
                                                    const char                      * address_string,
                                                    int                               len_address_string);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_load_mirror_(struct hmpprti_grouplet * const & grouplet,
                                                    const char                      * address_string,
                                                    const int                       & async,
                                                    const int                       * slice_count,
                                                    ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_wait_load_mirror_(struct hmpprti_grouplet * const & grouplet,
                                                         const char                      * address_string,
                                                         int                               len_address_string);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_store_mirror_(struct hmpprti_grouplet * const & grouplet,
                                                     const char                      * address_string,
                                                     const int                       & async,
                                                     const int                       * slice_count,
                                                     ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_wait_store_mirror_(struct hmpprti_grouplet * const & grouplet,
                                                          const char                      * address_string,
                                                          int                               len_address_string);

HMPPRTI_API
void hmpprti_run_directive_(struct hmpprti_grouplet * const & grouplet,
                            const char                      * file_name,
                            const int                       & line_number,
                            const int                       * semantic,
                            ...);

HMPPRTI_API
void hmpprti_enqueue_directive_(struct hmpprti_grouplet * const & grouplet,
                                const char                      * file_name,
                                const int                       & line_number,
                                const int                       * semantic,
                                ...);

HMPPRTI_API
void hmpprti_run_queue_(void);

#ifdef __cplusplus
}
#endif

#endif
