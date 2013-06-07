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
  INTEGER, PARAMETER :: acc_device_kind = 4 ! Must be equal to C_INT
  INTEGER(acc_device_kind), PARAMETER :: acc_device_none = 0 !< no device
  INTEGER(acc_device_kind), PARAMETER :: acc_device_default = 1 !< default device type
  INTEGER(acc_device_kind), PARAMETER :: acc_device_host = 2 !< host device
  INTEGER(acc_device_kind), PARAMETER :: acc_device_not_host = 3 !< not host device
  INTEGER(acc_device_kind), PARAMETER :: acc_device_cuda = 4 !< CUDA device
  INTEGER(acc_device_kind), PARAMETER :: acc_device_opencl = 5 !< OpenCL device
  INTERFACE
    FUNCTION acc_get_num_devices( devicetype ) result(out)
      USE ISO_C_BINDING
      INTEGER(C_INT) :: out
      INTEGER(C_INT), INTENT(IN) :: devicetype
    END FUNCTION acc_get_num_devices
    SUBROUTINE acc_set_device_type( devicetype )
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: devicetype
    END SUBROUTINE acc_set_device_type
    FUNCTION acc_get_device_type() result(out)
      USE ISO_C_BINDING
      INTEGER(C_INT) :: out
    END FUNCTION acc_get_device_type
    SUBROUTINE acc_set_device_num( devicenum, devicetype )
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: devicenum
      INTEGER(C_INT), INTENT(IN) :: devicetype
    END SUBROUTINE acc_set_device_num
    FUNCTION acc_get_device_num( devicetype ) result(out)
      USE ISO_C_BINDING
      INTEGER(C_INT) :: out
      INTEGER(C_INT), INTENT(IN) :: devicetype
    END FUNCTION acc_get_device_num
    ! Async
    FUNCTION acc_async_test( arg ) result(out)
      USE ISO_C_BINDING
      LOGICAL(C_BOOL) :: out
      INTEGER(C_INT), INTENT(IN) :: arg
    END FUNCTION acc_async_test
    FUNCTION acc_async_test_all() result(out)
      USE ISO_C_BINDING
      LOGICAL(C_BOOL) :: out
    END FUNCTION acc_async_test_all
    SUBROUTINE acc_async_wait( arg )
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: arg
    END SUBROUTINE acc_async_wait
    SUBROUTINE acc_async_wait_all()
    END SUBROUTINE acc_async_wait_all
    !> init
    SUBROUTINE acc_init( devicetype )
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: devicetype
    END SUBROUTINE acc_init
    SUBROUTINE acc_shutdown( devicetype )
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: devicetype
    END SUBROUTINE acc_shutdown
    FUNCTION acc_on_device( devicetype ) result(out)
      USE ISO_C_BINDING
      LOGICAL(C_BOOL) :: out
      INTEGER(C_INT), INTENT(IN) :: devicetype
    END FUNCTION acc_on_device
  END INTERFACE
