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
/// \internal

#include <stdint.h>

namespace hmpprt
{

/// \internal
struct cipher_t;

#ifndef _WIN32
typedef ::uint8_t  u08;
typedef ::uint16_t u16;
typedef ::uint32_t u32;
typedef ::uint64_t u64;

typedef ::int8_t   s08;
typedef ::int16_t  s16;
typedef ::int32_t  s32;
typedef ::int64_t  s64;
#else
typedef unsigned __int8   u08;
typedef unsigned __int16  u16;
typedef unsigned __int32  u32;
typedef unsigned __int64  u64;

typedef signed __int8   s08;
typedef signed __int16  s16;
typedef signed __int32  s32;
typedef signed __int64  s64;
#endif

} // namespace hmpprt
