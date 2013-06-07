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

#ifndef __INTEL_COMPILER
#error "Intel compiler not detected"
#endif

#if __INTEL_COMPILER < 1000
#error "Intel compiler < 10.00 are not supported"
#endif

#define KIND_INT8   C_INT8_T
#define KIND_INT16  C_INT16_T
#define KIND_INT32  C_INT32_T
#define KIND_INT64  C_INT64_T

#define KIND_REAL32 C_FLOAT
#define KIND_REAL64 C_DOUBLE

#define BOOL_TYPE LOGICAL(4)

! All ifort compilers use -1 to represent .TRUE.
#define TRUE_IS_MINUS_1

! CLASS(*) is supported since ifort 11.01
!#if __INTEL_COMPILER >= 1101
!#define CLASS_SUPPORT
!#endif

#define DEC_NO_ARG_CHECK_SUPPORT

#include <hmpprt/hmpprt_module.F90>

