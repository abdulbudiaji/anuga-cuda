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
#ifndef HMPPABI_SIGNATURE_H
#define HMPPABI_SIGNATURE_H

#include <hmppabi/Type.h>

#include <string>
#include <vector>
#include <map>

namespace hmppabi
{

class Signature
{
public:
  Signature();
  Signature(const Signature & signature);
  ~Signature();

  void setFunctionName(const std::string & name);
  std::string getFunctionName() { return function_name_; }

  std::string getParameterName(int index) const;
  void addParameter(const Type & type);
  void addParameter(const std::string & name, const Type & type);
  Type & getParameter(const std::string & name);
  Type & getParameter(int index) { return parameters_[index]; }
  const Type & getParameter(int index) const { return parameters_[index]; }
  MemorySpace getParameterPointedMemorySpace(int index);
  MemorySpace getParameterPointedMemorySpace(const std::string & name);
  Intent getParameterPointedIntent(int index);
  int getParameterIndex(const std::string & name);
  int getNumberOfParameters() const;

  Signature bind(const std::map<std::string, int> & symbols);

  bool hasReturnType();
  Type & getReturnType();
  void setReturnType(const Type & type);

  static Signature from_dpil_text(const std::string & s);
  static Signature from_c_text(const std::string & s);
  static Signature from_fortran_text(const std::string & s);
  std::string asDPILText();
  std::string asCText();

private:
  std::vector<Type> parameters_;
  Type return_type_;
  std::string function_name_;
};

} // namespace hmppabi

#endif
