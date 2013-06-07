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
#ifndef OPENACCI_DATA_H
#define OPENACCI_DATA_H

#include <openacci/Mirror.h>
#include <openacci/DeviceResident.h>

#include <hmpprt/hmpprt.h>
#include <cassert>

namespace openacci
{


/// Represent a piece of data allocated on the device ;
/// Kind of decorator of hmpprt::Data.
struct Data
{
  Data(void * host_address, size_t size);

  Data(void           * host_address,
       size_t           start,
       size_t           length,
       size_t           elemsize);
  ~Data();

  void useBufferOf(Data * parent);
  void useSubBufferOf(Data * parent);
  void useDevicePtr(hmpprt::Device * device, hmpprt::MemorySpace, void * addr);
  void updateDeviceResidentShape(DeviceResident dr);
  void useDeviceResident(DeviceResident dr);

  hmpprt::Data * getBuffer() { return buffer_; }

  size_t getSize() { return buffer_->getSize(); }
  void * getDeviceAddress() { return buffer_->getAddress(); }
  size_t getOffset() { return buffer_->getOffset(); }
  hmpprt::MemorySpace getMemorySpace() { return buffer_->getMemorySpace(); }
  hmpprt::Device * getDevice() { return buffer_->getDevice(); }

  void allocate(const char          * file_name,
                int                   line_number,
                hmpprt::Device      * device,
                hmpprt::MemorySpace   memory_space,
                hmpprt::Queue       * queue,
                const char          * qdesc);

  void free(const char    * file_name,
            int             line_number,
            hmpprt::Queue * queue,
            const char    * qdesc);

  void dispose(hmpprt::Queue * queue);

  void upload(const char    * file_name,
              int             line_number,
              hmpprt::Queue * queue,
              const char    * qdesc);

  void download(const char    * file_name,
                int             line_number,
                hmpprt::Queue * queue,
                const char    * qdesc);

  void partialUpload(const char    * file_name,
              int             line_number,
                size_t offset,
                size_t size,
              hmpprt::Queue * queue,
              const char    * qdesc);

  void partialDownload(const char    * file_name,
                int             line_number,
                size_t offset,
                size_t size,
                hmpprt::Queue * queue,
                const char    * qdesc);

  bool contains(void * host_address, size_t size);

  void setVariableName(const char * s)
  {
    variable_name_ = s;
  }

  void setMirror(MirrorMap::iterator mi) { mi_ = mi; }
  MirrorMap::iterator getMirror() { return mi_; }

  // Reference host address, passed in argument to the region
  void * getHostAddress() { return host_address_; }

  // Beginning of the mirrorred data
  void * getStartHostAddress() { return (char *) host_address_ + getOffset(); }

  const char * getVariableName() { return variable_name_.c_str(); }
  size_t getNumElements() { return length_; }
  size_t getFirstElement() { return start_; }
  size_t getElementSize() { return element_size_; }

  bool           pushed;             //< Data to pop from the mirror stack when leaving the region
  bool           allocated;          //< Data to free when leaving region
  bool           copyin;             //< Data to copy to host when calling the region
  bool           copyout;            //< Data to copy back to host when leaving the region
  bool           device_resident;
  bool           value() { return buffer_ == 0; }

private:
  //@{
  /// For debugging purpose
  std::string    variable_name_;
  size_t         start_;
  size_t         length_;
  size_t         element_size_;
  //@}

  /// host_address_ + device_data_.getOffset() is mirrorred at device_data_.getAddress()
  void         * host_address_;      

  MirrorMap::iterator
                 mi_;                //< Mirror referencing this data

  hmpprt::Data * buffer_;
  Data * parent_data_;
};


} // openacci namespace

#endif // OPENACCI_DATA_H
