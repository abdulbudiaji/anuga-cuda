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

#ifndef HMPPRT_OPENCL_GRID_H
#define HMPPRT_OPENCL_GRID_H

#include <string>

#include <hmpprt/Common.h>
/* #include <hmpprt/OpenCLTypes.h> */
#include <hmpprt/OpenCLDevice.h>
#include <hmpprt/HostTypes.h>
#include <hmpprt/Context.h>

namespace hmpprt
{

class OpenCLModule;

/// \internal
class OpenCLGrid
{

public:

  HMPPRT_API
  OpenCLGrid(OpenCLModule * module, const char * function_name);

  HMPPRT_API
  ~OpenCLGrid();

  HMPPRT_API
  OpenCLModule * getModule() const;

  HMPPRT_API
  const std::string & getFunctionName() const;

private:

  OpenCLModule  * module_;
  std::string   function_name_;

};

class OpenCLGridCall
{
  friend class OpenCLDevice;

public:

  HMPPRT_API
  OpenCLGridCall();

  HMPPRT_API
  ~OpenCLGridCall();

public:

  // local parameters
  HMPPRT_API
  void addLocalParameter(MemorySpace ms, void * data, const char * name);
  HMPPRT_API
  void addLocalParameter(s08 value, const char * name);
  HMPPRT_API
  void addLocalParameter(s16 value, const char * name);
  HMPPRT_API
  void addLocalParameter(s32 value, const char * name);
  HMPPRT_API
  void addLocalParameter(float value, const char * name);
  HMPPRT_API
  void addLocalParameter(double value, const char * name);
  HMPPRT_API
  void addLocalParameter(const void * value_ptr, size_t value_size, const char * name);
  template <MemorySpace MS, typename BT>
  void addLocalParameter(DevicePtr<MS,BT> data, const char * name)
  {
    addLocalParameter(MS, data.addr, name);
  }

  // const paremeters
  HMPPRT_API
  void addConstParameter(const void * host_ptr, const char * name);
  HMPPRT_API
  void addConstParameter(MemorySpace ms, const void * host_ptr, const char * name);
  template <MemorySpace MS, typename BT>
  void addConstParameter(DevicePtr<MS,BT> dp, const char * name)
  {
    addConstParameter(MS, dp.addr, name);
  }

  // shared parameters
  HMPPRT_API
  void addSharedParameter(MemorySpace ms, void * data, const char * name);
  template <MemorySpace MS, typename BT>
  void addSharedParameter(DevicePtr<MS,BT> data, const char * name)
  {
    addSharedParameter(MS, data.addr, name);
  }

  HMPPRT_API
  void setWorkDim(int work_dim) { work_dim_ = work_dim; }
  HMPPRT_API
  void setBlockSizeX(int block_size_x) { block_size_x_ = block_size_x; }
  HMPPRT_API
  void setBlockSizeY(int block_size_y) { block_size_y_ = block_size_y; }
  HMPPRT_API
  void setBlockSizeZ(int block_size_z) { block_size_z_ = block_size_z; }
  HMPPRT_API
  void setSizeX(int grid_size_x) { grid_size_x_ = grid_size_x; }
  HMPPRT_API
  void setSizeY(int grid_size_y) { grid_size_y_ = grid_size_y; }
  HMPPRT_API
  void launch(OpenCLGrid * grid, OpenCLDevice * device);

//private:

  std::vector<void *> param_list_;
  std::map<std::string, void *> constant_pointers_;
  std::map<std::string, void *> constant_values_;

  enum ParamType { PT_DATA_PTR, PT_CHAR, PT_SHORT, PT_INT, PT_FLOAT, PT_DOUBLE, PT_SIZED_VALUE };
  std::vector<ParamType> param_types_;
  int shared_size_;
  int shared_reduction_size_;

  int work_dim_;
  int block_size_x_;
  int block_size_y_;
  int block_size_z_;
  int grid_size_x_;
  int grid_size_y_;

};

} // namespace hmpprt

#endif // HMPPRT_OPENCL_GRID_H
