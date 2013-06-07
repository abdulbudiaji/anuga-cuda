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
#ifndef _HMPPABI_PARSER_H_
#define _HMPPABI_PARSER_H_

#include <hmpperr/Error.h>

#include <hmppabi/Signature.h>

namespace hmppabi
{
  struct ParsingResult
  {
  public:
    ParsingResult()
      {}
    ParsingResult(Type type, int eaten_chars) :
      type_(type), eaten_chars_(eaten_chars) {}
    ParsingResult(Signature signature, int eaten_chars) :
      signature_(signature), eaten_chars_(eaten_chars) {}
    ~ParsingResult()
      {}

  public:
    Signature getSignature()
      { return signature_; }
    Type getType()
      { return type_; }
    int getEatenChars()
      { return eaten_chars_; }

  private:
    Type type_;
    Signature signature_;
    int eaten_chars_;
  };


  ParsingResult parse_dpil_signature(const std::string & signature);
  ParsingResult parse_dpil_signature(const std::string & s, const std::map<std::string, Type> & compounds);

  ParsingResult parse_dpil_type(const std::string & type);
  ParsingResult parse_dpil_type(const std::string & type, const std::map<std::string, Type> & compounds);


  ParsingResult parse_c_signature(const std::string & s);
  ParsingResult parse_c_signature(const std::string & s, const std::map<std::string, Type> & compounds);

  ParsingResult parse_c_type(const std::string & type);
  ParsingResult parse_c_type(const std::string & type, const std::map<std::string, Type> & compounds);


  ParsingResult parse_fortran_signature(const std::string & s);
  ParsingResult parse_fortran_signature(const std::string & s, const std::map<std::string, Type> & compounds);

  ParsingResult parse_fortran_type(const std::string & type);
  ParsingResult parse_fortran_type(const std::string & type, const std::map<std::string, Type> &);

  std::vector<Type> parse_fortran_declare(const std::string & s);

} // namespace hmppabi

#endif // _HMPPABI_PARSER_H_
