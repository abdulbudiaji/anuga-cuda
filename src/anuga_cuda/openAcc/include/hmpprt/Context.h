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
/// \internal

#ifndef HMPPRT_CONTEXT_H
#define HMPPRT_CONTEXT_H

#include <map>
#include <string>
#include <cstddef>
#include <iostream>

#include <hmpprt/Common.h>
#include <hmpprt/Data.h>
#include <hmpprt/Device.h>

namespace hmpprt
{

class CUDADevice;
class OpenCLDevice;
class Grouplet;


template <MemorySpace MS, typename T>
struct DeviceAccess;

template <MemorySpace MS, typename T>
struct DevicePtr;

template <MemorySpace MS, typename T>
struct DeviceAccess;

template <MemorySpace MS, typename T>
struct DevicePtr
{
  T * addr;

  explicit DevicePtr(T * value = 0)
    : addr(value)
  { }

  inline DevicePtr<MS,T> operator = (void * host_address);

  // Use of generated copy cosntructor and copy assign operator

  DevicePtr<MS,T> operator ++ ()
  {
    addr += 1;
    return * this;
  }
  DevicePtr<MS,T> operator -- ()
  {
    addr -= 1;
    return * this;
  }
  DevicePtr<MS,T> operator ++ (int)
  {
    DevicePtr<MS,T> p(*this);
    addr += 1;
    return p;
  }
  DevicePtr<MS,T> operator -- (int)
  {
    DevicePtr<MS,T> p(*this);
    addr -= 1;
    return p;
  }

#define TYPE(OFFSET_TYPE)                                                   \
  DevicePtr<MS,T> operator += (OFFSET_TYPE offset)                          \
  {                                                                         \
    addr += offset;                                                         \
    return * this;                                                          \
  }                                                                         \
  DevicePtr<MS,T> operator -= (OFFSET_TYPE offset)                          \
  {                                                                         \
    addr -= offset;                                                         \
    return * this;                                                          \
  }                                                                         \
  friend DevicePtr<MS,T> operator + (DevicePtr<MS,T> p, OFFSET_TYPE offset) \
  {                                                                         \
    p.addr += offset;                                                       \
    return p;                                                               \
  }                                                                         \
  friend DevicePtr<MS,T> operator - (DevicePtr<MS,T> p, OFFSET_TYPE offset) \
  {                                                                         \
    p.addr -= offset;                                                       \
    return p;                                                               \
  }                                                                         \
  friend DevicePtr<MS,T> operator + (OFFSET_TYPE offset, DevicePtr<MS,T> p) \
  {                                                                         \
    p.addr += offset;                                                       \
    return p;                                                               \
  }
  
  TYPE(signed char)
  TYPE(unsigned char)
  TYPE(short int)
  TYPE(int)
  TYPE(long int)
  TYPE(long long int)
  TYPE(unsigned short int)
  TYPE(unsigned int)
  TYPE(unsigned long int)
  TYPE(unsigned long long int)

#undef TYPE

  operator bool() 
  {
    return addr != 0 ; 
  }
  friend ptrdiff_t operator - (DevicePtr<MS,T> left, DevicePtr<MS,T> right)
  {
    return left.addr - right.addr;
  }
  friend bool operator == (DevicePtr<MS,T> left, DevicePtr<MS,T> right)
  {
    return left.addr == right.addr;
  }
  friend bool operator != (DevicePtr<MS,T> left, DevicePtr<MS,T> right)
  {
    return left.addr != right.addr;
  }

  DeviceAccess<MS,T> operator * ()
  {
    return operator[](0);
  }

  inline DeviceAccess<MS,T> operator [] (size_t index);

  inline DeviceAccess<MS,T> operator -> ();

  T * host() ; 
  
};

template <MemorySpace MS, typename T>
T * DevicePtr<MS,T>::host() 
{ 
  return  ((*this)[0]).read() ; 
} 

/// \internal
template <MemorySpace MS, typename T>
struct DeviceAccess
{
  DeviceAccess(T * addr, size_t index)
    : addr(addr), index(index)
  {
  }

  T * addr;
  size_t index;

  DevicePtr<MS,T> operator & ()
  {
    DevicePtr<MS,T> ptr;
    ptr.addr = addr + index;
    return ptr;
  }

  operator T ()
  {
    return * read();
  }

  DeviceAccess<MS,T> operator = (T t)
  {
    write(t);
    return * this;
  }

  inline void write(T t);
  inline T * read();

  DeviceAccess<MS,T> operator = (DeviceAccess<MS,T> other)
  {
    return operator=((T) other);
  }

  T * operator -> ()
  {
    return read();
  }
};


/// \internal
class Context : NonCopyable
{
public:
  HMPPRT_API
  static Context * getInstance();

  Context();
  ~Context();

public:
  HMPPRT_API
  void setGrouplet(Grouplet * grouplet);
  HMPPRT_API
  Grouplet * getGrouplet();
  HMPPRT_API
  void setDevice(Device * device);

  // Get current device
  // Only valid during/inside a codelet call
  HMPPRT_API
  Device * getDevice();

  // Get current device and check this is a CUDADevice
  // Only valid during/inside a codelet call
  HMPPRT_API
  CUDADevice * getCUDADevice();

  // Get current device and check this is a OpenCLDevice
  // Only valid during/inside a codelet call
  HMPPRT_API
  OpenCLDevice * getOpenCLDevice();

  // For use from Codelet::call (through ArgumentList)
  HMPPRT_API
  void addData(Data * data);

  HMPPRT_API
  void * read_from_host(Data *);

  //@{
  // Create/free a temporary on device
  HMPPRT_API
  void allocate(void **src, MemorySpace memory_space, size_t size);
  HMPPRT_API
  void free(void ** src);
  //@}

  // Upload and invalidate every temporary buffer on host side
  HMPPRT_API
  void invalidate();
  HMPPRT_API
  std::map<void*,Data*>::iterator getDataFromHostAddress(void * a);
  HMPPRT_API
  std::map<void*,Data*>::iterator getDataFromDeviceAddress(void * a);

  HMPPRT_API
  Queue * getQueue() { return queue_; }
  HMPPRT_API
  void setQueue(Queue * queue) { queue_ = queue; }

private:
  static __thread Context * instance_;

private:
  Grouplet * grouplet_;
  Device * device_;
  Queue * queue_;

  std::map<void*,Data*> host_map_;
  std::map<void*,Data*> device_map_;
};

// DevicePtr

template <MemorySpace MS, typename T>
inline DeviceAccess<MS,T> DevicePtr<MS,T>::operator [] (size_t index)
{
  return DeviceAccess<MS,T>(addr, index);
}

template <MemorySpace MS, typename T>
inline DevicePtr<MS,T> DevicePtr<MS,T>::operator = (void * host_address)
{
  Context * ctx = Context::getInstance();
  std::map<void*,Data*>::iterator it = ctx->getDataFromHostAddress(host_address);
  void * ha_base = it->first;
  size_t offset = ((char *) host_address) - ((char *) ha_base);
  Data * d = it->second;
  addr = (T *) (((char *) d->getAddress()) + offset);
  return * this;
}

template <MemorySpace MS, typename T>
DeviceAccess<MS,T> DevicePtr<MS,T>::operator -> ()
{
  return operator[](0);
}


// DeviceAccess

template <MemorySpace MS, typename T>
inline void DeviceAccess<MS,T>::write(T t)
{
  * read() = t;
}

template <MemorySpace MS, typename T>
inline T * DeviceAccess<MS,T>::read()
{
  Context * ctx = Context::getInstance();
  std::map<void*,Data*>::iterator it = ctx->getDataFromDeviceAddress(addr);
  void * addr_base = it->first;
  size_t offset = ((char *)addr) - ((char *)addr_base);
  Data * d = it->second;
  T * host_addr = (T *) ctx->read_from_host(d);

  return (T *) (((char *)(host_addr + index)) + offset);
}

} // namespace hmpprt

#endif // HMPPRT_CONTEXT_H
