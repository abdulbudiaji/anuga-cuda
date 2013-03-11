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
#ifndef OPENACCI_DEVICE_RESIDENT_H
#define OPENACCI_DEVICE_RESIDENT_H

#include <stdlib.h>

namespace openacci
{
  struct Data;

  class DeviceResident
  {
  public:
    DeviceResident(const void * address)
      : array_desc_addr_((size_t *) address)
    { }
  
    size_t   getNbDim() { return array_desc_addr_[1]; }

    size_t   getElemSize() { return array_desc_addr_[2]; }

    size_t * getBounds() { return &array_desc_addr_[3]; }

    size_t getNumElements();

    void setHandle(Data * data) { array_desc_addr_[17] = (size_t) data; }

    Data * getHandle() { return (Data *)array_desc_addr_[17]; }

    bool isAllocated() { return (array_desc_addr_[0] != 0); }
    void setAllocated(bool allocated = true)
    {
      // FIXME: I guess this is useless
      array_desc_addr_[0] = (size_t) allocated;
    }
  
  private:
    size_t * array_desc_addr_;
  };
}

#endif // OPENACCI_DEVICE_RESIDENT_H
