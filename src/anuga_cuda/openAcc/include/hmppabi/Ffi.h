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
#ifndef HMPPABI_FFI_H
#define HMPPABI_FFI_H

#include <hmppabi/Signature.h>
#include <hmppabi/Type.h>


namespace hmppabi
{

  class Ffi
  {
  public:
    Ffi(Signature * signature);
    void call(void (* function_ptr)(), int arguments_count, void * args[]);
    void * translateType(Type & type);
    ~Ffi();

  private:

    std::vector<void *> ffi_types_;
    void *              ffi_cif_;
    Signature *         signature_;

    std::vector<void *>  store_struct_type;
    std::vector<void **> store_struct_type_elements;
  };

} // namespace hmppabi

#endif
