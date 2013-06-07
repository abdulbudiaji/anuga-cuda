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
#ifndef OPENACCI_OPENACCI_H
#define OPENACCI_OPENACCI_H

#include <openacc/openacc.h>
#include <hmpprt/Target.h>

#ifdef OPENACCI_API
#  undef OPENACCI_API
#endif /* OPENACCI_API */

#ifdef _WIN32
#  ifdef OPENACCI_BUILD
#    define OPENACCI_API __declspec(dllexport)
#  else /* ! OPENACCI_BUILD */
#    ifndef OPENACC_BUILD
#      pragma comment(lib, "openacci")
#    endif
#    define OPENACCI_API __declspec(dllimport)
#  endif /* OPENACCI_BUILD */
#else /* ! _WIN32 */
#  define OPENACCI_API
#endif /* _WIN32 */

#include <openacci/Context.h>

namespace openacci
{

OPENACCI_API
const char * get_device_type_string(acc_device_t device_type);

OPENACCI_API
hmpprt::Target get_target(acc_device_t device_type);

enum AsyncMode
{
    SYNC       = 0,
    QUEUE      = 1,
    FULL_ASYNC = 2
};

}

#endif
