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
#ifndef HMPPRT_MALLOC_MANAGER_H_
#define HMPPRT_MALLOC_MANAGER_H_

#include <hmpprt/Device.h>
#include <hmpprt/Context.h>

#include <hmppmem/hmppmem.h>

#include <hmppos/hmppos.h>

namespace hmpprt
{

  struct MallocManager
  {
    HmppMemEnableAllocationHook enable_alloc_hook_func_;
    HmppMemEnableFreeHook       enable_free_hook_func_;
    HmppMemMallocFunc           mem_host_alloc_;
    HmppMemFreeFunc             mem_free_host_;

    HmppMemDisableHook          disable_hooks_func_;
    Device *                    device_;

    bool                        enabled_;

    MallocManager();

    static MallocManager * getInstance()
    {
      static MallocManager * instance_;
      hmppos::singleton(&instance_);
      return instance_;
    }

    void setDevice(Device * d) { device_ = d; }
    Device * getDevice()       { return device_; }

    void loadHmppmemLibrary();

    ~MallocManager()
    {
      disableHooks();
    }

    void enableHooks(HmppMemMallocFunc memHostAlloc, HmppMemFreeFunc memFreeHost);
    bool disableHooks();
  };


  namespace 
  {
    struct HandleHookDisabler
    {
      HandleHookDisabler()
      {
        alloc_hook_ = MallocManager::getInstance()->mem_host_alloc_;
        free_hook_ = MallocManager::getInstance()->mem_free_host_;
        malloc_was_enabled_ = MallocManager::getInstance()->disableHooks();
      }
      
      ~HandleHookDisabler()
      {
        if (malloc_was_enabled_)
        {
          MallocManager::getInstance()->enableHooks(*alloc_hook_, *free_hook_);
        }
      }

      bool malloc_was_enabled_;
      HmppMemMallocFunc alloc_hook_;
      HmppMemFreeFunc free_hook_;
    };
  }

} // namespace hmpprt

#endif
