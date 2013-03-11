!
! Copyright (C) 2008-2013 CAPS entreprise.  All Rights Reserved.
! 
! The source code contained or described herein and all documents  related
! to the source code ("Material") are owned  by  CAPS  entreprise  or  its
! suppliers or licensors.
! 
! Title to the Material remains with CAPS entreprise or its suppliers  and
! licensors.  The Material contains trade secrets and proprietary and con-
! fidential information of CAPS entreprise or its suppliers and licensors.
! 
! The Material is protected by the French intellectual property code,  in-
! tellectual property laws and international treaties.  No part of the Ma-
! terial may be used, copied, reproduced, modified,  published,  uploaded,
! posted, transmitted, distributed or disclosed in any  way  without  CAPS
! entreprise's prior express written permission.
! 
! No license under any patent, copyright, trade secret or other  intellec-
! tual property right is granted to or conferred upon you by disclosure or
! delivery of the Material, either expressly, by implication,  inducement,
! estoppel or otherwise.
! 
! Any license under such intellectual property rights  must  be  expressed
! and approved by CAPS entreprise in writing.
!

#ifndef __GNUC__
#error "GFortran compiler not detected"
#endif

#ifndef _LANGUAGE_FORTRAN
#error "GFortran compiler not detected"
#endif

#define VERSION(major,minor,patch)   (major*10000+minor*100+patch)

#define GNU_VERSION VERSION(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__)

! GFortran compilers older than 4.3 probably work fine but
! they lack a lot of features (e.g. ISO_C_BINDING) and
! they are also quite immature.
#if ( GNU_VERSION < VERSION(4,3,0) )
#error "GFortran compilers older than 4.3.0 are not supported"
#endif

#define KIND_INT8   C_INT8_T
#define KIND_INT16  C_INT16_T
#define KIND_INT32  C_INT32_T
#define KIND_INT64  C_INT64_T

#define KIND_REAL32 C_FLOAT
#define KIND_REAL64 C_DOUBLE

#define BOOL_TYPE LOGICAL(4)

! All gfortran compilers use +1 to represent .TRUE.
#define TRUE_IS_PLUS_1

! CLASS(*) is still not fully supported by gfortran 4.6
!#if ( GNU_VERSION >= VERSION(4,7,0) )
!#define CLASS_SUPPORT
!#endif

#undef VERSION
#undef GNU_VERSION

#include "hmpprt_module.F90"

