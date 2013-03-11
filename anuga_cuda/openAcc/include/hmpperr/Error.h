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
#ifndef HMPPERR_ERROR_H
#define HMPPERR_ERROR_H

#include <string>
#include <exception>
#include <vector>

// Remark: The known errors are defined in Error.def


namespace hmpperr
{

/// Base class for every exception arisen by HMPP RT
class Error : public std::exception
{

public:

  Error(const char * file, int line) throw ();
  Error() throw ();

  virtual ~Error() throw ();

  virtual Error * clone() const throw () __attribute__ ((warn_unused_result)) = 0;

  virtual void rethrow() const = 0;

  virtual const char * getClassName() const throw () = 0;

  virtual const char * what() const throw () __attribute__ ((warn_unused_result));

  virtual int id() const = 0;

  void setText(const char * fmt, ...) throw () __attribute__((format(printf, 2, 3)));
  void setText() throw ();

private:

  std::string what_;

};

#define HMPP_ERROR_DEF(klass,Kode,Id)                 \
class klass : public Error                            \
{                                                     \
public:                                               \
  klass(const char * file, int line) throw ();        \
  klass() throw ();                                   \
  virtual ~klass() throw ();                          \
  virtual klass * clone() const throw ();             \
  virtual void rethrow() const;                       \
  virtual const char * getClassName() const throw (); \
  virtual int id() const ;                            \
};
#include <hmpperr/Error.def>
#undef HMPP_ERROR_DEF


#if defined(RELEASE) || defined(NDEBUG)

#define HMPPERR_THROW(klass, ...) \
do                                \
{                                 \
  hmpperr::klass error___;        \
  error___.setText(__VA_ARGS__);  \
  throw error___;                 \
}                                 \
while (false)

#else

#define HMPPERR_THROW(klass, ...)              \
do                                             \
{                                              \
  hmpperr::klass error___(__FILE__, __LINE__); \
  error___.setText(__VA_ARGS__);               \
  throw error___;                              \
}                                              \
while (false)

#endif

std::vector<std::string> get_error_names();


}

#endif
