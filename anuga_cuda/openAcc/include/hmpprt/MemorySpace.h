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
#ifndef HMPPRT_MEMORY_SPACE_H
#define HMPPRT_MEMORY_SPACE_H

#include <hmpprt/Common.h>

namespace hmpprt
{

/// Memory space
enum MemorySpace
{
#define D(id, en, sn) en = id,
#include <hmpprt/MemorySpace.def>
#undef D
  MS_MAX
};

/// \return Memory space name
HMPPRT_API
const char * get_memory_space_name(MemorySpace memory_space);

} // namespace hmpprt

#endif // HMPPRT_MEMORY_SPACE_H