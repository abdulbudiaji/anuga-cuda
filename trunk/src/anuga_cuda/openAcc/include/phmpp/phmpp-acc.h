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
#ifndef PHMPP_ACC_H
#define PHMPP_ACC_H

#include <phmpp/phmpp.h>

#ifdef __cplusplus
extern "C"
{
#endif //__cplusplus

#define PHMPP_ACC_VERSION           100
#define PHMPP_ACC_BASE            30000
#define PHMPP_ACC_MASK             0xFF

  /* Gives the version of the profling semantic extension. */
  static inline int phmpp_get_acc_version(void)           { return PHMPP_ACC_VERSION; }

  /* Gives the base of the profling semantic extension. */
  static inline int phmpp_get_acc_base(void)              { return PHMPP_ACC_BASE; }

  /* Gives the mask of the profling semantic extension. */
  static inline int phmpp_get_acc_mask(void)              { return PHMPP_ACC_MASK; }

  /* Return a value different from 0 if the given semnatic is a phmpp_acc_event_t */
  static inline int phmpp_is_acc(unsigned int semantic)   { return ((semantic | PHMPP_ACC_MASK) == (PHMPP_ACC_BASE | PHMPP_ACC_MASK)); }
  
  typedef enum  {
    PHMPP_ACC_UNDEFINED = PHMPP_ACC_BASE,
#define D(NAME, EVENT, RTI, DESC) \
    PHMPP_ACC_ ## NAME,
#include <phmpp/acc.def>
#undef D
    PHMPP_ACC_NOMORE
  } phmpp_acc_event_t;

  PHMPP_API
  const char * phmpp_get_acc_event_name(unsigned int semantic);
  PHMPP_API
  const char * phmpp_get_acc_event_description(unsigned int semantic);

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //PHMPP_ACC_H
