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
#ifndef HMPPRT_GROUPLET_H
#define HMPPRT_GROUPLET_H

#include <string>
#include <vector>
#include <map>

#include <hmpprt/Common.h>
#include <hmpprt/DeviceManager.h>

namespace hmpprt
{

class Codelet;

/// List of codelets
typedef std::vector<Codelet *> CodeletList;

/// A grouplet is a set of Codelet functions generated for the same target.
/// \see Codelet
class Grouplet : private NonCopyable
{

public:

  /// Load a grouplet
  /// \param file_name is the path of the dynamic loadable library which contains the generated code of codelets.
  HMPPRT_API
  explicit Grouplet(const std::string & file_name);

  /// Load a grouplet
  /// \param source_file is the path of a source file (which contains at least one hmpp entrypoint pragma.)
  /// \param target is the target for what the source file has been generated with HMPP.
  HMPPRT_API
  Grouplet(const std::string & source_file, Target target);

  HMPPRT_API
  ~Grouplet();

  /// \return the file name of the grouplet.
  HMPPRT_API
  const std::string & getFileName() const;

  /// \return the target of the grouplet.
  HMPPRT_API
  Target getTarget() const;

  /// \return the list of codelets present in this grouplet.
  HMPPRT_API
  CodeletList listCodelets();

  /// \return the Codelet belonging to this Grouplet with the given function name.
  HMPPRT_API
  Codelet * getCodeletByName(const std::string & codelet_name);

public:

  /// \internal
  HMPPRT_API
  void setTarget(Target target);

  /// \internal
  HMPPRT_API
  void addSignature(const char * codelet_name, const char * signature_string);

  /// \internal
  HMPPRT_API
  bool hasSignature(const std::string & codelet_name);

private:

  typedef std::map<std::string, Codelet *> CodeletsByName;
  typedef std::vector<void *> SignatureList;

  SignatureList signatures_;

  std::string      file_name_;
  void           * handle_;
  Target           target_;
  CodeletList      codelets_;
  CodeletsByName   codelets_by_name_;

  void load();
  void unload();

};

} // namespace hmpprt

#endif // HMPPRT_GROUPLET_H
