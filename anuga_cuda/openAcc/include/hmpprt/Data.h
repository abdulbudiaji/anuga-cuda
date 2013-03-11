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
#ifndef HMPPRT_DATA_H
#define HMPPRT_DATA_H

#include <hmpprt/Common.h>
#include <hmpprt/NonCopyable.h>
#include <hmpprt/MemorySpace.h>
#include <cstring>

namespace hmpprt
{

class Device;
class Queue;

/// A Data object encapsulate a buffer of memory that resides either in host memory
/// or in the memory of a Device.
class Data
{
public:
  /// Create a sparse data object residing on the given device.
  /// \param device is the device on which the data can be allocated
  /// \param address is the device address
  /// \param size is in bytes
  /// \param offset is subtracted to the device address at codelet call
  /// \param memory_space is used to choose the type of memory to use on the device.
  /// \param host_address is the host address of the data.
  HMPPRT_API
  Data(Device * device, void * address, size_t size, size_t offset, MemorySpace memory_space, void * host_address=0);

  /// Create a sparse data object residing on the given device.
  /// \param device is the device on which the data can be allocated
  /// \param size is in bytes
  /// \param offset is subtracted to the device address at codelet call
  /// \param memory_space is used to choose the type of memory to use on the device.
  /// \param host_address is the host address of the data.
  HMPPRT_API
  Data(Device * device, size_t size, size_t offset, MemorySpace memory_space, void * host_address=0);

  /// Create a data object residing on the given device.
  /// \param device is the device on which the data can be allocated
  /// \param size is in bytes
  /// \param memory_space is used to choose the type of memory to use on the device.
  /// \param host_address is the host address of the data.
  HMPPRT_API
  Data(Device * device, size_t size, MemorySpace memory_space, void * host_address=0);

  /// Create a data object residing on the given device.
  /// \param device is the device on which the data can be allocated
  /// \param size is in bytes
  /// The default memory space of the given device is used.
  /// \param host_address is the host address of the data.
  HMPPRT_API
  Data(Device * device, size_t size, size_t offset = 0, void * host_address=0);

  HMPPRT_API
  Data(const Data &);

  HMPPRT_API
  ~Data();

  /// \return the device associated with this data.
  HMPPRT_API
  Device * getDevice();

  /// \internal
  HMPPRT_API
  void setDevice(Device * device) { device_ = device; }

  /// \return the size of this data, in bytes.
  HMPPRT_API
  size_t getSize() const;

  /// Set a new size for this data, in bytes.
  HMPPRT_API
  void setSize(size_t size) { size_ = size; }

  /// \return the memory space of this data.
  HMPPRT_API
  MemorySpace getMemorySpace() const;

  /// Set a new memory space name
  HMPPRT_API
  void setMemorySpace(MemorySpace ms) { memory_space_ = ms; }

  /// \return the name of the memory space of this data (for debugging purpose).
  HMPPRT_API
  const char * getMemorySpaceName() const;

  // Get the device address of the data.
  // For advance usage only.
  /// \internal
  HMPPRT_API
  void * getAddress() const;

  // Set the device address of the data.
  // For advance usage only.
  /// \internal
  HMPPRT_API
  void setAddress(void * device_address);

  /// Allocates the data in the device memory.
  HMPPRT_API
  void allocate();

  /// Frees the allocated device memory.
  HMPPRT_API
  void free();

  /// Uploads data from given host address.
  HMPPRT_API
  void upload(const void * host_address, Queue * queue = 0);

  /// Partial upload of data from the given host address.
  /// Transfers 'size' bytes from 'host_address' + 'offset', and write them at 'offset'.
  /// Other bytes keep their previous value.
  HMPPRT_API
  void upload(const void * host_address, size_t offset, size_t size, Queue * queue = 0);

  /// Downloads data into the given host address.
  HMPPRT_API
  void download(void * host_address, Queue * queue = 0);

  /// Partial dowload of data into the given host address.
  /// Transfers 'size' bytes from 'offset', and write them at 'host_address' + 'offset'.
  /// Other bytes keep their previous value.
  HMPPRT_API
  void download(void * host_address, size_t offset, size_t size, Queue * queue = 0);

  // Invalidated by a call to a codelet taking this data as argument.
  //@{
  /// \internal
  HMPPRT_API
  void * getHandle() { return handle_; }
  HMPPRT_API
  void setHandle(void * h) { handle_ = h; }
  //@}

  //@{
  /// \internal
  HMPPRT_API
  void * getHostAddress() { return host_address_; }
  HMPPRT_API
  void setHostAddress(void * host_address) { host_address_ = host_address; }
  //@}

  // \return the offset to subtract to device_address_ before call
  // to compensate the accesses in the targetted code.
  /// \internal
  HMPPRT_API
  size_t getOffset() { return offset_; }

  /// \internal
  HMPPRT_API
  void setOffset(size_t offset) { offset_ = offset; }

protected:

  Device      * device_;
  size_t        size_;
  size_t        offset_;
  MemorySpace   memory_space_;
  void        * device_address_;
  void        * host_address_;
  void        * handle_;

};

} // namespace hmpprt

#endif // HMPPRT_DATA_H
