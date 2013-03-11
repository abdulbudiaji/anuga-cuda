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
#ifndef HMPPRT_HOST_DEVICE_H
#define HMPPRT_HOST_DEVICE_H

#include <hmpprt/Device.h>

namespace hmpprt
{

/// The HostDevice represents the access to the main central processing units.
/// \see Device
class HostDevice : public Device
{
public:
  /// \internal
  static void create_devices(DeviceList &);

  virtual ~HostDevice();

  virtual MemorySpace getDefaultMemorySpace() const;

  virtual bool hasDoublePrecision() const;

  virtual bool hasHostMemoryIntegration() const;

  /// This methods are not implemented for the host device.
  virtual unsigned int getDeviceMaxClockRate() const;

  void allocate(Data * data);

  void free(Data * data);

  void upload(Data * data, const void * host_address, size_t offset, size_t size, Queue * q = 0);

  void download(Data * data, void * host_address, size_t offset, size_t size, Queue * q = 0);

  virtual void releaseHandle(void * handle);

  virtual void synchronize(Queue * queue = 0);

  virtual bool isPageLocked(void * host_address, size_t size);

  virtual void replaceAllocs(int substitute);

  bool checkDeviceData(void * device_address, void * host_address, size_t size, Queue * queue);

  /// \internal
  virtual void createContext();

protected:
  HostDevice(unsigned long       id,
             const std::string & label);
};

} // namespace hmpprt

#endif // HMPPRT_HOST_DEVICE_H
