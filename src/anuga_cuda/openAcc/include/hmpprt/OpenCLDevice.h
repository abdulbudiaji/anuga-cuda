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
#ifndef HMPPRT_OPENCL_DEVICE_H
#define HMPPRT_OPENCL_DEVICE_H

#include <map>

#include <hmpprt/HostDevice.h>

namespace hmpprt
{

class OpenCLGrid;

class OpenCLGridCall;

/// OpenCL-capable device
/// \see Device
class OpenCLDevice : public HostDevice
{
  friend class OpenCLContextLocker;

public:

  /// \internal
  static void initialize();

  /// \internal
  static void create_devices(DeviceList &);

  virtual ~OpenCLDevice();

  virtual MemorySpace getDefaultMemorySpace() const;

  virtual bool hasDoublePrecision() const;

  virtual bool hasHostMemoryIntegration() const;

  virtual unsigned int getDeviceMaxClockRate() const;

  //@{
  /// Access to information specific to OpenCL devices.
  const std::string & getOpenCLDeviceName() const;

  const std::string & getOpenCLDeviceVendor() const;

  const std::string & getOpenCLDriverVersion() const;

  const std::string & getOpenCLDeviceVersion() const;

  size_t getOpenCLDeviceGlobalMemorySize() const;
  //@}

  void allocate(Data * data);

  void free(Data * data);

  void upload(Data * data, const void * host_address, size_t offset, size_t size, Queue * q = 0);

  void download(Data * data, void * host_address, size_t offset, size_t size, Queue * q = 0);

  //@{
  /// \internal
  void launchGrid(OpenCLGrid * grid, OpenCLGridCall * call, Queue * queue);
  void releaseHandle(void * handle);
  void * getHandle();
  void * getTransferQueue(Queue * queue);
  void * getKernelQueue(Queue * queue);
  void replaceAllocs(int);
  void * getZeroCopyBuffer(void * mapped);
  void addZeroCopyBuffer(void * mapped, void * buffer);

  /// Performs a synchronous download and compare the downloaded data 
  /// to the data at address host_address
  bool checkDeviceData(void * device_address, void * host_address, size_t size, Queue * queue);
  //@}

  /// \internal
  virtual void createContext();

private:

  void * handle_;

  OpenCLDevice(unsigned long       id,
               const std::string & label,
               void              * handle);

  void constructBuildOptions();

  std::map<void *, void*> zerocopy_map_;

};

/// Lock the OpenCL context of a device for use from the calling thread.
/// For advance use (ie. interoperabiliity of hmpprt with hand-written OpenCL code)
/// \see OpenCLDevice
class OpenCLContextLocker
{

public:

  /// Lock the context of a device for exclusive use from the current thread.
  /// \param device is the device on which we want the OpenCL code to run.
  explicit OpenCLContextLocker(OpenCLDevice * device);

  /// Unlock the context of the device
  ~OpenCLContextLocker();

private:

  OpenCLDevice * device_;

};

} // namespace hmpprt

#endif // HMPPRT_OPENCL_DEVICE_H
