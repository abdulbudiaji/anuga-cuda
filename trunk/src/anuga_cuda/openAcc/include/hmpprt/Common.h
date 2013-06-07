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
#ifndef HMPPRT_COMMON_H
#define HMPPRT_COMMON_H

#ifdef HMPPRT_API
#  undef HMPPRT_API
#endif /* HMPPRT_API */

#ifdef _WIN32

#  ifdef HMPPRT_BUILD
#    define HMPPRT_API __declspec(dllexport)
#  else /* ! HMPPRT_BUILD */
#    pragma comment(lib, "hmpprt")
#    define HMPPRT_API __declspec(dllimport)
#  endif /* HMPPRT_BUILD */

#  define __attribute__(ARG)
#  define __thread __declspec(thread)
#  include <BaseTsd.h>
#  define ssize_t SSIZE_T

#else /* ! _WIN32 */

#  define HMPPRT_API

#endif /* _WIN32 */


#include <string>
#include <vector>

namespace hmpplog
{
  class Logger;
}

namespace hmpprt
{

#define HMPPRT_ENVIRONMENT_DEF(VAR, TYPE, DEFAULT) \
extern const TYPE VAR;
#include <hmpprt/environment.def>
#undef HMPPRT_ENVIRONMENT_DEF

  struct LoggerWrapper
  {
    hmpplog::Logger * logger_;

    LoggerWrapper();

    ~LoggerWrapper();
  };

  void setup_logger(LoggerWrapper ** w);

  HMPPRT_API
  inline hmpplog::Logger * logger()
  {
    static LoggerWrapper * instance_ = 0;

    if (!instance_)
    {
      setup_logger(&instance_);
    }

    return instance_->logger_;
  }

class AutoIndent
{
public:
  HMPPRT_API
  AutoIndent(int indent = 2)
    : inc_(indent)
  {
    indent_ += inc_;
  }
  
  HMPPRT_API
  ~AutoIndent()
  {
    indent_ -= inc_;
  }

  HMPPRT_API
  static void inc(int indent = 2)
  {
    indent_ += indent;
  }

  HMPPRT_API
  static void dec(int indent = 2)
  {
    indent_ -= indent;
  }

  HMPPRT_API
  static int get()
  {
    return indent_;
  }

private:
  static __thread int indent_;

  int inc_;
};

}

#define HMPPRT_DEBUG(fmt, ...)  do { HMPPLOG_LOG(hmpprt::logger(), DEBUG, "%*s" fmt, hmpprt::AutoIndent::get(), "" , ## __VA_ARGS__); } while (0)
#define HMPPRT_INFO(fmt, ...)   do { HMPPLOG_LOG(hmpprt::logger(), INFO , "%*s" fmt, hmpprt::AutoIndent::get(), "" , ## __VA_ARGS__); } while (0)
#define HMPPRT_WARN(fmt, ...)   do { HMPPLOG_LOG(hmpprt::logger(), WARN , "%*s" fmt, hmpprt::AutoIndent::get(), "" , ## __VA_ARGS__); } while (0)
#define HMPPRT_ERROR(fmt, ...)  do { HMPPLOG_LOG(hmpprt::logger(), ERROR, "%*s" fmt, hmpprt::AutoIndent::get(), "" , ## __VA_ARGS__); } while (0)
#define HMPPRT_FATAL(fmt, ...)  do { HMPPLOG_LOG(hmpprt::logger(), FATAL, "%*s" fmt, hmpprt::AutoIndent::get(), "" , ## __VA_ARGS__); } while (0)

#ifdef _WIN32
HMPPRT_API
void setCrashHandler();
#  define HMPPRT_THROW(err, fmt, ...) do { setCrashHandler(); HMPPRT_ERROR(fmt, __VA_ARGS__); HMPPERR_THROW(err, fmt, __VA_ARGS__); } while(0)
#else
#  define HMPPRT_THROW(err, ...) HMPPERR_THROW(err, __VA_ARGS__)
#endif

#endif /* HMPPRT_COMMON_H */
