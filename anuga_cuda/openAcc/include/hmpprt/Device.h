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
#ifndef HMPPRT_DEVICE_H
#define HMPPRT_DEVICE_H

#include <string>
#include <vector>

#include <hmpprt/Common.h>
#include <hmpprt/Target.h>
#include <hmpprt/NonCopyable.h>
#include <hmpprt/MemorySpace.h>

namespace hmpprt
{

class HostDevice;
class CUDADevice;
class OpenCLDevice;
class Codelet;
class ArgumentList;
class Device;
class Queue;
class Data;

/// Return a readable name for the given target
HMPPRT_API
const char * get_target_name(Target target);

HMPPRT_API
Target get_target_from_string(const std::string & target);

// raise an error or print a warning depending on environment 
// variable HMPPRT_CHECK_TRANSFERS's value
void raise_transfer_error_warning(void * src, void * dst, bool uploading);

/// List of Device objects
typedef std::vector<hmpprt::Device *> DeviceList;

/// Base class for a Hardware Accelerator Device which can runs codelets.
class Device : private NonCopyable
{
public:

  HMPPRT_API
  virtual ~Device();

  /// Acquire this device for exclusive use by the calling thread.
  /// \throw hmpperr::LockFailError if this device is already in use.
  /// \throw hmpperr::IOError if an error occured when communicating with the resource manager.
  /// \see release()
  HMPPRT_API
  void acquire();

  /// Try to acquire the Device.
  /// \return true if succeeds.
  HMPPRT_API
  bool tryAcquire();

  /// Release this device.
  /// \see acquire()
  HMPPRT_API
  void release();

  virtual void createContext() = 0;

  /// Execute the given codelet on this Device, using the given list of arguments
  HMPPRT_API
  void call(Codelet * codelet, ArgumentList & arguments, Queue * q = 0);

  /// \return the default memoryspace for this device.
  HMPPRT_API
  virtual MemorySpace getDefaultMemorySpace() const = 0;

  /// \return trus if double floating point precision is supported by this device.
  HMPPRT_API
  virtual bool hasDoublePrecision() const = 0;

  /// \return true if the device memory is shared with the host.
  HMPPRT_API
  virtual bool hasHostMemoryIntegration() const = 0;

  /// \return the max clock rate of the device.
  HMPPRT_API
  virtual unsigned int getDeviceMaxClockRate() const = 0;

  /// \return the unique ID of the device.
  HMPPRT_API
  unsigned long getId() const;

  /// \return a human readable description of this device
  HMPPRT_API
  const std::string & getLabel() const;

  /// Allocate memory of the given memory space.
  /// \param data data to allocate
  HMPPRT_API
  virtual void allocate(Data * data) = 0;

  /// Free memory of the given memory space.
  /// \param data data to free
  HMPPRT_API
  virtual void free(Data * data) = 0;

  /// Upload data from host to device.
  /// \param data data to transfer
  /// \param host_address the host address of data to transfer
  /// \param offset the offset to start the transfer
  /// \param size the size of data to transfer
  /// \param q queue where operation is performed
  HMPPRT_API
  virtual void upload(Data * data, const void * host_address, size_t offset, size_t size, Queue * q = 0) = 0;

  /// Download data from host to device.
  /// \param data data to transfer
  /// \param host_address the host address of data to transfer
  /// \param offset the offset to start the transfer
  /// \param size the size of data to transfer
  /// \param q queue where operation is performed  HMPPRT_API
  HMPPRT_API
  virtual void download(Data * data, void * host_address, size_t offset, size_t size, Queue * q = 0) = 0;

  /// Copy between data objects
  HMPPRT_API
  virtual void copy(Data * dst_data, Data * src_data, Queue * q = 0);

  /// \internal
  virtual void releaseHandle(void * handle) = 0;

  /// \internal
  virtual void synchronize(Queue * queue = 0) = 0;

  /// \internal
  virtual bool isPageLocked(void * host_address, size_t size) = 0;

  /// \internal
  virtual void replaceAllocs(int) = 0;

  /// Check whether the device data is equal to the host data
  /// returns true if the device_data is equal to data located at host_address
  /// else, returns false
  virtual bool checkDeviceData(void * device_address, void * host_address, size_t size, Queue * queue) = 0;

protected:
  Device(unsigned long       id,
         const std::string & label);

private:
  unsigned long id_;
  std::string   label_;
};

} // namespace hmpprt

#endif // HMPPRT_DEVICE_H
