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
#ifndef HMPPRT_ARGUMENT_LIST_H
#define HMPPRT_ARGUMENT_LIST_H

#include <vector>
#include <list>
#include <map>
#include <string>

#include <hmpprt/Common.h>
#include <hmpprt/MemorySpace.h>

namespace hmpprt
{

class Codelet;
class Data;
class Context;

/// This class provides a temporary object used to bind Data objects and scalar
/// values as a list of parameters to call a codelet.
/// Use generated default and copy constructor and copy operator.
class ArgumentList
{
public:
  HMPPRT_API
  ArgumentList();
  HMPPRT_API
  ArgumentList(const ArgumentList &);
  HMPPRT_API
  ~ArgumentList();

  //@{
  /// Enqueue the given Data object in the argument list.
  HMPPRT_API
  void addDataArgument(hmpprt::Data * data) { addArgument().data = data; }
  HMPPRT_API
  void addDataArgument(hmpprt::Data & data) { addDataArgument(& data); }
  HMPPRT_API
  void addArgument(hmpprt::Data * data) { addDataArgument(data); }
  HMPPRT_API
  void addArgument(hmpprt::Data & data) { addDataArgument(&data); }
  //@}

  //@{
  /// Add a scalar passed by value to the argument list.
  HMPPRT_API
  void addArgumentByCopy(const void*addr,size_t size);
  template <typename T>
  void addArgument(const T & n) { addArgumentByCopy(&n, sizeof(T)); }
  //@}

  /// \return the number of arguments
  HMPPRT_API
  int getNumberOfArguments() const { return args_count_; }

  /// \internal
  HMPPRT_API
  void getArgumentValues(void * args[]);

  /// \internal
  HMPPRT_API
  void mapData(Context * context);

  /// \internal
  HMPPRT_API
  Data * getDataAt(int i) { return getArgumentAt(i).data; }

  /// \internal
  HMPPRT_API
  void addArgumentFrom(const ArgumentList & source, int i) { source.getArgumentAt(i).clone(addArgument()); }

private:
  struct Argument
  {
    void clone(Argument &) const;

    Data * data;
    char * value;
    size_t size;
  };

  const Argument & getArgumentAt(int i) const { return ptArgs()[i]; }
  Argument & addArgument();

  static const int ARGS_COUNT = 10;
  Argument st_args_[ARGS_COUNT];
  int args_count_;
  std::vector<Argument> dn_args_;

  const Argument * ptArgs() const { return args_count_ <= ARGS_COUNT ? st_args_ : dn_args_.data(); }
  Argument       * ptArgs()       { return args_count_ <= ARGS_COUNT ? st_args_ : dn_args_.data(); }

  ArgumentList & operator = (const ArgumentList &);
};


} // namespace hmpprt

#endif // HMPPRT_ARGUMENT_LIST
