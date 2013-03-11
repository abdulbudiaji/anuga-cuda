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
#ifndef HMPPRT_CODELET_H
#define HMPPRT_CODELET_H

#include <string>

#include <hmpprt/Common.h>
#include <hmpprt/NonCopyable.h>
#include <hmpprt/MemorySpace.h>
#include <hmpprt/Intent.h>
#include <hmpprt/Queue.h>

namespace hmpprt
{

class Grouplet;
class Device;
class ArgumentList;

/// A Codelet is pure function which has been generated to run on a Device.
/// Several Codelet objects coming from the same source file compiled for a
/// specific Target are grouped into a Grouplet.
/// \see Device
/// \see Grouplet
/// \see Target
class Codelet : private NonCopyable
{
  friend class Grouplet;

  HMPPRT_API
  Codelet(Grouplet * grouplet, const std::string & name, void (* function)(), void * signature);

public:
  HMPPRT_API
  ~Codelet();

  /// \return the Grouplet this Codelet belongs to.
  HMPPRT_API
  Grouplet * getGrouplet();

  /// Return the memoryspace of the argument at given index.
  /// Useful to allocate a Data object suitable for the call.
  /// \param index starts at zero for the first argument.
  HMPPRT_API
  MemorySpace getMemorySpaceByIndex(int index);

  /// Return the memoryspace of the argument of given name.
  /// Useful to allocate a Data object suitable for the call.
  /// \param name is the formal name of the argument.
  HMPPRT_API
  MemorySpace getMemorySpaceByName(const std::string & name);

  /// Return the intent of the argument at given index.
  /// \param index starts at zero for the first argument.
  HMPPRT_API
  Intent getIntentByIndex(int index);

  /// \return the function name of this Codelet.
  HMPPRT_API
  const std::string & getName() const;

  /// Invoke the codelet with the given device and arguments
  HMPPRT_API
  void call(Device * device, ArgumentList & arguments, Queue * queue = 0);

  /// \return the number of parameters to this codelet
  HMPPRT_API
  int getNumberOfParameters() const;

  /// \param size returns the size of a parameter
  /// \return true if the parameter is passed by value to this codelet
  HMPPRT_API
  bool getValueSizeByIndex(int index, size_t & size) const;

  /// \return the name of a parameter
  HMPPRT_API
  std::string getNameByIndex(int index) const;

private:
  Grouplet     * grouplet_;
  std::string    name_;
  void        (* function_)();
  void         * signature_;
  void         * ffi_;

};

} // namespace hmpprt

#endif // HMPPRT_CODELET_H
