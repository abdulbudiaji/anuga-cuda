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
#ifndef HMPPRT_DEVICE_MANAGER_H
#define HMPPRT_DEVICE_MANAGER_H

#include <hmpprt/Common.h>
#include <hmpprt/Device.h>
#include <map>

namespace hmpprt
{


/// Class DeviceManager allows to retrieve available devices on the system.
/// This class is a singleton.
class DeviceManager : private NonCopyable
{
public:
  /// \return the DeviceManager singleton
  HMPPRT_API
  static DeviceManager * getInstance();

  /// Destroy all the created devices and cleanup the allocated memory
  HMPPRT_API
  void cleanUp();

  /// \internal
  static bool hasInstance();

  /// Hardware introspection.
  /// \param output receives the list of devices suitable for the specified target.
  /// \param target can be ANY_TARGET to get devices of any target.
  HMPPRT_API
  void searchDevices(Target target, DeviceList & output);

  /// Acquire for exclusive use the first 'number' available devices which supports the given target.
  /// \param output receives the acquired devices
  /// \param number if the number of device to acquire
  /// \param target can be ANY_TARGET to get devices of any target.
  /// \throw DeviceError if not enough suitable devices are available
  HMPPRT_API
  void acquireDevices(Target target, size_t number, DeviceList & output);

  /// \return the host device singleton
  HMPPRT_API
  Device * getHostDevice();

  /// Acquire for exclusive use the first available CUDA device.
  /// \throw DeviceError if not enough suitable devices are available
  HMPPRT_API
  Device * acquireCUDADevice();

  /// Acquire for exclusive use the first available CUDA device.
  /// \throw DeviceError if not enough suitable devices are available
  HMPPRT_API
  Device * acquireOpenCLDevice();

  /// Acquire for exclusive use the first available device for the given target
  HMPPRT_API
  Device * acquireDevice(Target target);

private:

  /// \internal
  std::map<Target, DeviceList> devices_;

  /// \internal
  void * mutex_;

public: // internal

  DeviceManager();

  ~DeviceManager();

};

} // namespace hmpprt

#endif // HMPPRT_DEVICE_MANAGER_H
