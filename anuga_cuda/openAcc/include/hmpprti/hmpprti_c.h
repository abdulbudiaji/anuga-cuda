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
#ifndef HMPPRTI_HMPPRTI_C_H
#define HMPPRTI_HMPPRTI_C_H

#include <stddef.h>

#ifdef HMPPRTI_API
#  undef HMPPRTI_API
#endif /* HMPPRTI_API */

#ifdef _WIN32
#  ifdef HMPPRTI_BUILD
#    define HMPPRTI_API __declspec(dllexport)
#  else /* ! HMPPRTI_BUILD */
#    pragma comment(lib, "hmpprti")
#    define HMPPRTI_API __declspec(dllimport)
#  endif /* HMPPRTI_BUILD */
#else /* ! _WIN32 */
#  define HMPPRTI_API
#endif /* _WIN32 */

#ifdef __cplusplus
extern "C"
{
#endif

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_lookup_grouplet(const char * grouplet_name);

HMPPRTI_API
void hmpprti_set_property(const char * name,
                          const char * value);

HMPPRTI_API
int hmpprti_is_initialized(const struct hmpprti_grouplet * grouplet);

HMPPRTI_API
void hmpprti_set_initialized(struct hmpprti_grouplet * grouplet,
                             int                       initialized);

HMPPRTI_API
void hmpprti_append_target(struct hmpprti_grouplet * grouplet,
                           const char              * target_name);

HMPPRTI_API
void hmpprti_append_buffer(struct hmpprti_grouplet * grouplet,
                           const char              * buffer_name,
                           size_t                    element_size,
                           int                       dimension_count,
                           int                       dynamic,
                           int                       lazy_transfer_at_callsite,
                           int                       load_once_between_two_allocates);

HMPPRTI_API
void hmpprti_append_codelet(struct hmpprti_grouplet * grouplet,
                            const char              * codelet_name,
                            const char              * function_name);

HMPPRTI_API
void hmpprti_append_buffer_to_codelet(struct hmpprti_grouplet * grouplet,
                                      int                       codelet_index,
                                      int                       buffer_index,
                                      int                       load_before_callsite,
                                      int                       store_after_callsite);

HMPPRTI_API
void hmpprti_append_buffer_to_codelet_hmpp2(struct hmpprti_grouplet * grouplet,
                                            int                       codelet_index,
                                            int                       buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_acquire(struct hmpprti_grouplet * grouplet,
                                               int                       default_device,
                                               int                       device,
                                               int                       exclusive);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_release(struct hmpprti_grouplet * grouplet,
                                               int                       default_device,
                                               int                       device);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_free_all(struct hmpprti_grouplet * grouplet);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_allocate(struct hmpprti_grouplet * grouplet,
                                                int                       buffer_index,
                                                int                       ignore_unbound,
                                                int                       default_device,
                                                int                       device,
                                                int                       dimension_count,
                                                ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_free(struct hmpprti_grouplet * grouplet,
                                            int                       buffer_index,
                                            int                       ignore_unbound);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_bind(struct hmpprti_grouplet * grouplet,
                                            int                       buffer_index,
                                            const char              * address_string,
                                            void                    * host_address);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_load(struct hmpprti_grouplet * grouplet,
                                            int                       buffer_index,
                                            int                       async,
                                            int                       slice_count,
                                            ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_wait_load(struct hmpprti_grouplet * grouplet,
                                                 int                       buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_store(struct hmpprti_grouplet * grouplet,
                                             int                       buffer_index,
                                             int                       async,
                                             int                       slice_count,
                                             ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_wait_store(struct hmpprti_grouplet * grouplet,
                                                  int                       buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_callsite(struct hmpprti_grouplet * grouplet,
                                                int                       codelet_index,
                                                int                       async);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_callsite_hmpp2(struct hmpprti_grouplet * grouplet,
                                                      int                       codelet_index,
                                                      int                       async,
                                                      int                       argument_count,
                                                      ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_standalone_callsite(struct hmpprti_grouplet * grouplet,
                                                           const char              * function_name,
                                                           int                       argument_count);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_synchronize(struct hmpprti_grouplet * grouplet,
                                                   int                       codelet_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_host_read(struct hmpprti_grouplet * grouplet,
                                                 int                       buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_host_write(struct hmpprti_grouplet * grouplet,
                                                  int                       buffer_index);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_allocate_mirror(struct hmpprti_grouplet * grouplet,
                                                       const char              * address_string,
                                                       int                       default_device,
                                                       int                       device,
                                                       size_t                    element_size,
                                                       int                       dimension_count,
                                                       ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_free_mirror(struct hmpprti_grouplet * grouplet,
                                                   const char              * address_string);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_load_mirror(struct hmpprti_grouplet * grouplet,
                                                   const char              * address_string,
                                                   int                       async,
                                                   int                       slice_count,
                                                   ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_wait_load_mirror(struct hmpprti_grouplet * grouplet,
                                                        const char              * address_string);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_store_mirror(struct hmpprti_grouplet * grouplet,
                                                    const char              * address_string,
                                                    int                       async,
                                                    int                       slice_count,
                                                    ...);

HMPPRTI_API
struct hmpprti_grouplet * hmpprti_push_wait_store_mirror(struct hmpprti_grouplet * grouplet,
                                                         const char              * address_string);

HMPPRTI_API
void hmpprti_run_directive(struct hmpprti_grouplet * grouplet,
                           int                     * error,
                           const char              * file_name,
                           int                       line_number,
                           int                       semantic,
                           ...);

HMPPRTI_API
void hmpprti_enqueue_directive(struct hmpprti_grouplet * grouplet,
                               const char              * file_name,
                               int                       line_number,
                               int                       semantic,
                               ...);

HMPPRTI_API
void hmpprti_run_queue(int * error);

HMPPRTI_API
void hmpprti_use_cuda_uva(int use);

#ifdef __cplusplus
}
#endif

#endif
