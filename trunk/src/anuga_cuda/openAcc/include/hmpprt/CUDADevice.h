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
#ifndef HMPPRT_CUDA_DEVICE_H
#define HMPPRT_CUDA_DEVICE_H

#include <map>

#include <hmpprt/HostDevice.h>
#include <hmpprt/HostTypes.h>

namespace hmpprt
{

class CUDAGrid;

class CUDAGridCall;

/// CUDA-capable device
/// \see Device
class CUDADevice : public HostDevice
{
  friend class CUDAContextLocker;

public:

  /// \internal
  static void initialize();

  /// \internal
  static void create_devices(DeviceList &);

  /// \internal
  virtual void createContext();

  virtual ~CUDADevice();

  virtual MemorySpace getDefaultMemorySpace() const;

  virtual bool hasDoublePrecision() const;

  virtual bool hasHostMemoryIntegration() const;

  virtual unsigned int getDeviceMaxClockRate() const;

  /// \return true is memory error correction code is enabled for this device.
  virtual bool hasECCEnabled() const;

  //@{
  /// Information specific for CUDA devices
  const std::string & getCUDADeviceName() const;

  int getCUDADriverVersion() const;

  void getCUDADeviceComputeCapability(int &major, int &minor) const;

  size_t getCUDADeviceTotalMemory() const;

  /// Return true if the CUDA context of this device is reusable for use from CUDA-runtime API.
  bool hasCUDAInteropCapabilities();
  //@}

  void allocate(Data * data);

  void free(Data * data);

  void upload(Data * data, const void * host_address, size_t offset, size_t size, Queue * q = 0);

  void download(Data * data, void * host_address, size_t offset, size_t size, Queue * q = 0);

  //@{
  /// \internal
  void launchGrid(CUDAGrid * grid, CUDAGridCall * call, Queue * q = 0);
  static void unloadModule(void * module);
  void lockContext();
  void unlockContext();
  void releaseHandle(void * handle);
  void synchronize(Queue * q = 0);
  void synchronizeWithoutSettingContext(Queue * q = 0);
  bool isPageLocked(void * host_address, size_t size);
  void replaceAllocs(int);

  /// Performs a synchronous download and compare the downloaded data 
  /// to the data at address host_address
  bool checkDeviceData(void * device_address, void * host_address, size_t size, Queue * queue);
  //@}


private:

  void uploadConstants(void * module, std::string function_name, std::map<std::string, void *> * pointers, std::map<std::string, void *> * values, Queue * q = 0);

  void * handle_;

  CUDADevice(unsigned long       id,
             const std::string & label,
             void              * handle);

};


/// Lock the CUDA context of a device for use from the calling thread.
/// For advance use (ie. interoperabiliity of hmpprt with hand-written CUDA code)
/// \see CUDADevice
class CUDAContextLocker
{

public:

  /// Lock the context of a device for exclusive use from the current thread.
  /// \param device is the device on which we want the CUDA code to run.
  explicit CUDAContextLocker(CUDADevice * device);

  /// Unlock the context of the device
  ~CUDAContextLocker();

  void unlock();

private:

  CUDADevice * device_;
  bool locked_;

};

} // namespace hmpprt

#endif // HMPPRT_CUDA_DEVICE_H
