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
#ifndef OPENACCI_REGION_H
#define OPENACCI_REGION_H

#include <openacci/Mirror.h>
#include <hmpprt/Queue.h>

namespace openacci
{
  class Context;

  /// Represents an acc region of code (data, kernels, parallel or the scope of a declare)
  struct Region
  {
    const char * getName() const;
    bool isKernels()  const;
    bool isParallel() const;
    bool isDeclare()  const;
    bool isData()     const;
    bool isHostData() const;
    
    int             kind;
    int             num_args;
    hmpprt::Queue * queue;
    std::string     qdesc;
    unsigned int    instance_id; //phmmp purpose

    // Add these information because we have no location if the region is global
    // (i.e declare in a module). 
    // In that case, the end of the region is not explicit, so has no location
    std::string     region_file_name;
    int             region_line_number;

    std::vector<Data> datas;
  };

} // openacci namespace

#endif // OPENACCI_REGION_H
