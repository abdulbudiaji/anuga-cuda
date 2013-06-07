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
#ifndef OPENACCI_MIRROR_H
#define OPENACCI_MIRROR_H

#include <vector>
#include <map>
#include <cstddef>

namespace openacci
{
  struct Data;

  /// A mirror is the name given to a host buffer which can be allocated and
  /// transfered several times on/to GPU (CPUBuffer would be a better name.)
  class Mirror
  {
  public:
    Mirror(size_t size) : fullsize_(size) { }

    /// \return a suitable Data object to use this mirror passing the given host address.
    Data * lookupData(void * host_address, size_t size);
  
    /// \return the size of the mirror
    size_t getSize() { return fullsize_; }

    Data * topData() { return data_stack_.back(); }
  
    void pushData(Data * data) { data_stack_.push_back(data); }
  
    bool empty() { return data_stack_.empty(); }
  
    void popData() { data_stack_.pop_back(); }
  
  private:
  
    size_t      fullsize_;             // Original size
    std::vector<Data*>   data_stack_;
  };

  typedef std::map<void *, Mirror> MirrorMap;
} // openacci namespace

#endif // OPENACCI_MIRROR_H
