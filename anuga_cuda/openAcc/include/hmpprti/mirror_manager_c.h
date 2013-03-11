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
#ifndef HMPPRTI_MIRROR_MANAGER_C_H
#define HMPPRTI_MIRROR_MANAGER_C_H


#include <stddef.h>
#include <hmpprt/hmpprt_c.h>
#include <hmpprti/hmpprti_c.h>


#ifdef __cplusplus
extern "C"
{
#endif

/// Mirror querying.
/// \param host_address is the source address of the mirror (host-side).
/// \return the hmpprt_data_t object related to this host address, or null if no mirror is found.
/// Each thread has its own mirror map.
/// The caller is responsible for making sure that the mirror is up-to-date at the time of use.
HMPPRTI_API
hmpprt_data_t hmpprti_get_mirror(const void * host_address);

/// Mirror and offset querying
/// \param host_address is the source address belonging to the mirror (host-side).
/// \param offset returns the offset between source address of the mirror and the host_address parameter
/// \return the hmpprt_data_t object related to this host address, or null if no mirror is found.
/// Each thread has its own mirror map.
/// The caller is responsible for making sure that the mirror is up-to-date at the time of use.
HMPPRTI_API
hmpprt_data_t hmpprti_get_mirror_and_offset(const void * host_address, size_t * offset);

/// Mirror registering
/// \param host_address is the source address associated with hmpprt_data object
/// \param data hmpprt_data object of the mirror
/// \return 1 if registering works, or 0 if mirror associated with host_address already exists
/// The caller is responsible for making sure that the hmpprt_data object is allocated.

HMPPRTI_API
int hmpprti_put_mirror(const void    * host_address,
                       hmpprt_data_t   data);

#ifdef __cplusplus
}
#endif

#endif

