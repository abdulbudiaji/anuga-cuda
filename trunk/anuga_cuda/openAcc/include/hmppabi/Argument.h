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
#ifndef HMPPABI_ARGUMENT_H
#define HMPPABI_ARGUMENT_H

#include <stdint.h>
#include <stdlib.h>

#include <hmppabi/Type.h>

namespace hmppabi
{

class Argument
{
public:
  Argument();
  Argument(void * content, Type type);
  Argument(const Argument & arg);
  ~Argument();

public:
#define D(num, name, ctype, dtyte, h5type, ffi_type, signed_, unsigned_) static Argument from_ ## dtype(ctype value);
#include <hmppabi/ScalarType.def>
#undef D

public:
  Type getType();
  void * getContent();
  void * getContentReference();

private:
  void * content_;
  Type type_;
};


} // namespace hmppabi

#endif



