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
#ifndef HMPPRTI_MIRROR_MANAGER_H
#define HMPPRTI_MIRROR_MANAGER_H

#include <string>

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

namespace hmpprt
{
class Data;
}

namespace hmpprti
{

HMPPRTI_API
bool try_get_mirror(const void     * host_address,
                    hmpprt::Data * & output);

HMPPRTI_API
bool try_get_mirror_and_offset(const void     * host_address,
                               hmpprt::Data * & output,
                               size_t         * offset);

HMPPRTI_API
bool try_put_mirror(const void   * host_address,
                    hmpprt::Data * data);

}

#endif
