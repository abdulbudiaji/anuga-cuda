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
!
! The following macros should be provided for each fortran compiler
!
!  KIND_INT8   = the INTEGER kind for 8bit integer
!  KIND_INT16  = the INTEGER kind for 16bit integer
!  KIND_INT32  = the INTEGER kind for 32bit integer
!  KIND_INT64  = the INTEGER kind for 64bit integer
!
!  KIND_REAL32 = the INTEGER kind for 32bit real
!  KIND_REAL64 = the INTEGER kind for 64bit real
!
!  TRUE_IS_MINUS_1 = define it if .TRUE. is -1
!  TRUE_IS_PLUS_1  = define it if .TRUE. is +1
!
!    ==> The .TRUE. value should be carrefully configured
!        to avoid strange errors
!
!    ==> Compile and execute the following code to
!        find out the .TRUE. value used by your compiler
!
!          PROGRAM test
!            LOGICAL(4) :: a = .TRUE.
!            CALL foo(a)
!          END PROGRAM test
!
!          SUBROUTINE foo(a)
!            INTEGER(4) :: a
!            PRINT *,".TRUE. = ",a
!          END SUBROUTINE foo
!
!  BOOL = the 32 bit LOGICAL type (e.g. LOGICAL(4) )
!
!  DEC_NO_ARG_CHECK_SUPPORT = only defined if the compiler supports
!                             the DEC directive NO_ARG_CHECK to disable
!                             type checking on some arguments
!
!  NO_ISO_C_BINDING = if defined then the module ISO_C_BINDING
!                     will not be used and a few missing macros
!                     such as C_SIZE_T, C_INT and C_INTPTR_T
!                     must be provided
!
!  DOC_MODE = should only be set to generate the Doxygen documentation
!


#ifndef DEC_NO_ARG_CHECK_SUPPORT
#warning "No support for DEC NO_ARG_CHECK directives. Some functions will be declared as EXTERNAL"
#endif



!> \defgroup grp_errors     Error Codes
!> @{
!> Below are all the possible values for the INTEGER(HMPPRT_ERROR_T) type used to
!> represent error status
!> @}
!> \defgroup grp_memspace   Memory Spaces
!> @{
!> Below are all the possible values for the INTEGER(HMPPRT_MEMSPACE_T) type used to
!> represent memory spaces
!> @}
!> \defgroup grp_integer    Integer Types
!> @{
!> The PARAMETERs specified below describe some INTEGER kind used
!> by the HMPPRT Fortran api
!> @}
!> \defgroup grp_targets    Targets
!> @{
!> Below are all the possible values for the INTEGER(HMPPRT_TARGET_T) type used to
!> represent the known HMPP targets
!> @}
!> \defgroup grp_handlers   Object Handlers
!> @{
!> Each of the compound types listed below represents an object manipulated
!> by the HMPP runtime
!>
!> The meaning of the field \c handle is not documented (this is actually
!> an address to the real object) but it can be assumed that a value of 0
!> represents an uninitialized object (e.g. NULL in C or C++)
!>
!> @}
!> \defgroup grp_device     Devices
!> @{
!> This section presents all functions related to device management
!> @}
!> \defgroup grp_grouplet   Grouplets
!> @{
!> This section presents all functions related to grouplet management
!> @}
!> \defgroup grp_codelet    Codelets
!> @{
!> This section presents all functions related to Codelet management
!> @}
!> \defgroup grp_arglist    Argument Lists
!> @{
!> This section presents all functions related to Argument List management
!> @}
!> \defgroup grp_queues     Execution Queues
!> @{
!> This section presents all functions related to Queue management
!> @}
!> \defgroup grp_data       Data Buffers
!> @{
!> This section presents all functions related to Data Buffer management
!> @}
!> \defgroup grp_cuda       CUDA Specific Functions
!> @{
!> @}
!> \defgroup grp_opencl     OpenCL Specific Functions
!> @{
!> @}
!>

!> \mainpage HMPPRT Fortran API
!>
!> \section intro_sec Introduction
!>
!> The HMPPRT Fortran API is provided via the module \ref hmpprt.
!>
!> Custom types and symbolic constants are provided by the module \ref hmpprt_enum which
!> should not be directly referenced since it is automatically imported by \ref hmpprt.
!>
!> Quick links:
!>   - \ref hmpprt      "The HMPPRT module"
!>   - \ref hmpprt_enum "The HMPPRT_ENUM module"
!>   - \ref grp_errors
!>   - \ref grp_memspace
!>   - \ref grp_integer
!>   - \ref grp_targets
!>   - \ref grp_handlers
!>   - \ref grp_device
!>   - \ref grp_grouplet
!>   - \ref grp_codelet
!>   - \ref grp_arglist
!>   - \ref grp_queues
!>   - \ref grp_data
!>   - \ref grp_cuda
!>   - \ref grp_opencl
!>
!>




!
! DOC_MODE is used to generate a better documentation in DOXYGEN
!
!
!
#ifdef DOC_MODE
#define DUMMY_INTERFACE
#define END_DUMMY_INTERFACE
#define NAMED_INTERFACE(x)
#define END_NAMED_INTERFACE(x)

! Provide a reference to name in another module
!
! For example: use REF(hmpprt_enum,HMPPRT_DEVICE_T) to refer to
!
#define REF(mod,name) \ref mod::name #name

! Convert x into "x"
#define TOSTRING(x)  #x

! Transform a function call into a reference
!
! The second argument should specify all arguments including the parenthesis be
!
! For example, REFCALL(xxx,(a,b,c)) will be transformed into a reference to xxx
! with the text "xxx(a,b,c)"
!
#define REFCALL(fun,args) \ref fun TOSTRING(fun##args)


#else
#define DUMMY_INTERFACE        INTERFACE
#define END_DUMMY_INTERFACE    END INTERFACE
#define NAMED_INTERFACE(x)     INTERFACE x
#define END_NAMED_INTERFACE(x) END INTERFACE x
#endif

!> \brief Internal module providing a few constants and types
MODULE hmpprt_enum

#ifndef NO_ISO_C_BINDING
  USE iso_c_binding
#endif

  IMPLICIT NONE


  !> Describe a HMPP Device
  !> \ingroup grp_handlers
  !> \ingroup grp_device
  TYPE HMPPRT_DEVICE_T
     INTEGER(C_INTPTR_T) :: handle !< Internal handle
  END TYPE HMPPRT_DEVICE_T

  !> Describe a HMPP Grouplet
  !> \ingroup grp_handlers
  !> \ingroup grp_grouplet
  TYPE HMPPRT_GROUPLET_T
     INTEGER(C_INTPTR_T) :: handle !< Internal handle
  END TYPE HMPPRT_GROUPLET_T

  !> Describe a HMPP Codelet
  !> \ingroup grp_handlers
  !> \ingroup grp_codelet
  TYPE HMPPRT_CODELET_T
     INTEGER(C_INTPTR_T) :: handle !< Internal handle
  END TYPE HMPPRT_CODELET_T

  !> Describe a HMPP Data buffer
  !> \ingroup grp_handlers
  !> \ingroup grp_data
  TYPE HMPPRT_DATA_T
     INTEGER(C_INTPTR_T) :: handle !< Internal handle
  END TYPE HMPPRT_DATA_T

  !> Describe a HMPP Argument list
  !> \ingroup grp_handlers
  !> \ingroup grp_arglist
  TYPE HMPPRT_ARGLIST_T
     INTEGER(C_INTPTR_T) :: handle !< Internal handle
  END TYPE HMPPRT_ARGLIST_T

  !> Describe a HMPP execution Queue
  !> \ingroup grp_handlers
  !> \ingroup grp_queue
  TYPE HMPPRT_QUEUE_T
     INTEGER(C_INTPTR_T) :: handle !< Internal handle
  END TYPE HMPPRT_QUEUE_T

  !> @}


  !> \addtogroup grp_integer
  !> @{

  !> An integer type used to represent any kind of address on a device.
  !>
  !> Its interpretation is device and memory space dependant. For example
  !>    - for a CUDA memory space, the integer value represents an address on the device
  !>    - for an OPENCL memory space, the integer value represents a cl_mem object
  !>    - for the HOST memory space, the integer value represents an address in the host main memory
  !>
  INTEGER, PARAMETER :: HMPPRT_ADDRESS_T = C_INTPTR_T

  !> The default integer type (should be 32 bit)
  INTEGER, PARAMETER :: HMPPRT_INT_T      = C_INT

  !> The integer type used to represent the size of an object (in bytes).
  !> This is typically an integer the size of a pointer
  INTEGER, PARAMETER :: HMPPRT_SIZE_T     = C_SIZE_T

  !> @}

  !> \addtogroup grp_targets
  !> @{

  !> The type INTEGER(HMPPRT_TARGET_T) describes the known HMPP targets
  INTEGER, PARAMETER :: HMPPRT_TARGET_T   = C_INT

  INTEGER(HMPPRT_TARGET_T), PARAMETER :: HMPPRT_ANY_TARGET    = 0   !< Not a real target. Can be used to select all targets
  INTEGER(HMPPRT_TARGET_T), PARAMETER :: HMPPRT_HOST_TARGET   = 1   !< The Host target (i.e. the CPU)
  INTEGER(HMPPRT_TARGET_T), PARAMETER :: HMPPRT_CUDA_TARGET   = 2   !< The CUDA target
  INTEGER(HMPPRT_TARGET_T), PARAMETER :: HMPPRT_OPENCL_TARGET = 4   !< The OPENCL target

  !> @}




  !> \addtogroup grp_memspace
  !> @{

  !> The type INTEGER(HMPPRT_MEMSPACE_T) describes the known HMPP memory spaces
  !>
  INTEGER, PARAMETER :: HMPPRT_MEMSPACE_T   = C_INT

  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_HOST           =    1   !< ...
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_CUDA_GLOB      = 1001   !< ...
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_CUDA_SHARED    = 1002   !< ...
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_CUDA_CONST     = 1003   !< ...
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_CUDA_PINNED    = 1004   !< ...
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_OPENCL_GLOB    = 2001   !< ...
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_OPENCL_SHARED  = 2002   !< ...
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_OPENCL_CONST   = 2003   !< ...
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_GHOST          =  601   !< internal
  INTEGER(HMPPRT_MEMSPACE_T), PARAMETER :: HMPPRT_MS_CIPHER         =  602   !< internal

  !> @}


  !> \addtogroup grp_errors
  !> @{

  !> The type INTEGER(HMPPRT_ERROR_T) describes the HMPP error codes
  INTEGER, PARAMETER :: HMPPRT_ERROR_T    = C_INT

  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_SUCCESS               =    0  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_FATAL_ERROR           =    1  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_SYSTEM_ERROR          =    2  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_IO_ERROR              =    3  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_LIBRARY_ERROR         =    4  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_LOOKUP_ERROR          =    5  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_NOT_IMPLEMENTED_ERROR =    6  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_VALUE_ERROR           =    7  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_DEVICE_ERROR          =    8  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_INDEX_ERROR           =    9  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_KEY_ERROR             =   10  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_DATASET_ERROR         =   11  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_SERIALIZATION_ERROR   =   12  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_PARSER_ERROR          =   13  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_CONTEXT_ERROR         =   14  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_TYPE_ERROR            =   15  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_LOCKFAIL_ERROR        =   16  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_ASSERTION_ERROR       =   17  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_UNKNOWN_ERROR         = 1001  !< ...
  INTEGER(HMPPRT_ERROR_T), PARAMETER :: HMPPRT_TIMEOUT_ERROR         =   -1  !< ...


  !> @}


#ifdef KIND_INT32
  !>
  !> A special value used as marker to help debug faulty formats
  !> (see hmpprt_arglist_fadd)
  !>
  INTEGER(KIND_INT32), PARAMETER :: HMPPRT_ARG_CHECK = 12345678
#endif

END MODULE hmpprt_enum

!> \brief This section provides all declarations needed to use the HMPPRT Fortran binding
MODULE hmpprt

  USE hmpprt_enum
  IMPLICIT NONE


#ifdef DOC_MODE

  !> \brief Generic subroutine to add arguments of many types
  !> \ingroup grp_arglist
  INTERFACE hmpprt_arglist_add_scalar
     MODULE PROCEDURE hmpprt_arglist_add_int8
     MODULE PROCEDURE hmpprt_arglist_add_int16
     MODULE PROCEDURE hmpprt_arglist_add_int32
     MODULE PROCEDURE hmpprt_arglist_add_int64
     MODULE PROCEDURE hmpprt_arglist_add_real32
     MODULE PROCEDURE hmpprt_arglist_add_complex32
     MODULE PROCEDURE hmpprt_arglist_add_real64
     MODULE PROCEDURE hmpprt_arglist_add_complex64
  END INTERFACE

  CONTAINS

#endif

  !>
  !> REMARK:
  !>
  !> Some functions and subroutines cannot be given a proper interface because
  !>    - they accept arguments of any types (i.e. the type 'void*' of C)
  !>        => The CLASS(*) feature of F2008 could work but it is not
  !>           supported by all compilers
  !>
  !>    - they accept a non-fixed number of arguments (i.e. the so-called 'varargs' of C)
  !>
  !> For each of those functions, an EXTERNAL declaration is given
  !> with a fake interface to describe the expected arguments
  !>

  ! ================================================================
  ! ==  Generic functions
  ! ================================================================

  DUMMY_INTERFACE
     !>
     !> Provide a readable name for a HMPP target (e.g. "HOST", "CUDA", ...)
     !>
     !> An empty name (only spaces) is returned if an unknown target is used.
     !>
     SUBROUTINE hmpprt_get_target_name(target,name)
       USE hmpprt_enum
       INTEGER(HMPPRT_TARGET_T) ,  INTENT(IN)  :: target  !< ...
       CHARACTER(*)             ,  INTENT(OUT) :: name !< ...
     END SUBROUTINE hmpprt_get_target_name
  END_DUMMY_INTERFACE

  ! ================================================================
  ! ==  HMPPRT_GROUPLET_T
  ! ================================================================

  DUMMY_INTERFACE
     !>
     !> Load a dynamic library generated by HMPP for a grouplet (or a single codelet)
     !>
     !> The filename should include the file extension (e.g. ".so" on Linux or ".dll"
     !> on windows)
     !>
     !> \ingroup grp_grouplet
     FUNCTION hmpprt_grouplet_load(grouplet,filename) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_GROUPLET_T), INTENT(OUT) :: grouplet  !< ...
       CHARACTER(*)              , INTENT(IN)  :: filename  !< ...
     END FUNCTION  hmpprt_grouplet_load
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Similar to hmpprt_grouplet_load but for a grouplet of a specific target.
     !>
     !> \ingroup grp_grouplet
     FUNCTION hmpprt_grouplet_load_t(grouplet,filename,target) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_GROUPLET_T), INTENT(OUT) :: grouplet  !< ...
       CHARACTER(*)              , INTENT(IN)  :: filename !< ...
       INTEGER(HMPPRT_TARGET_T) ,  INTENT(IN)  :: target  !< ...
     END FUNCTION hmpprt_grouplet_load_t
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Delete a grouplet previously allocated by hmpprt_grouplet_load or
     !> hmpprt_grouplet_load_t
     !>
     !> \ingroup grp_grouplet
     FUNCTION hmpprt_grouplet_delete(grouplet) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_GROUPLET_T), INTENT(IN)  :: grouplet  !< ...
     END FUNCTION  hmpprt_grouplet_delete
  END_DUMMY_INTERFACE

 DUMMY_INTERFACE
     !>
     !> Provide the filename associated to a grouplet
     !>
     !> \ingroup grp_grouplet
     SUBROUTINE hmpprt_grouplet_get_filename(grouplet,filename)
       USE hmpprt_enum
       TYPE(HMPPRT_GROUPLET_T), INTENT(IN)  :: grouplet  !< ...
       CHARACTER(*)              , INTENT(OUT) :: filename!< ...
     END SUBROUTINE hmpprt_grouplet_get_filename
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !> Provide the target for which a grouplet was compiled
     !> \ingroup grp_grouplet
     FUNCTION hmpprt_grouplet_get_target(grouplet) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_TARGET_T)                :: out
       TYPE(HMPPRT_GROUPLET_T), INTENT(IN)  :: grouplet  !< ...
     END FUNCTION hmpprt_grouplet_get_target
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !> Provide the number of codelets into a grouplet
     !> \ingroup grp_grouplet
     FUNCTION hmpprt_grouplet_get_codelet_count(grouplet) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                   :: out
       TYPE(HMPPRT_GROUPLET_T), INTENT(IN)  :: grouplet  !< ...
     END FUNCTION hmpprt_grouplet_get_codelet_count
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> Provide the codelet of the specified rank in a grouplet
     !> or HMPPRT_NO_CODELET if an invalid rank is specified
     !>
     !> Remark: The first codelet has rank 1
     !>
     !> \ingroup grp_grouplet
     FUNCTION hmpprt_grouplet_get_codelet(grouplet,codelet,rank) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_CODELET_T) , INTENT(OUT) :: codelet  !< ...
       TYPE(HMPPRT_GROUPLET_T), INTENT(IN)  :: grouplet  !< ...
       INTEGER(HMPPRT_INT_T)     , INTENT(IN)  :: rank  !< ...
     END FUNCTION hmpprt_grouplet_get_codelet
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Provide the codelet of the specified name in a grouplet
     !> or HMPPRT_NO_CODELET if an invalid name is specified
     !>
     !> \ingroup grp_grouplet
     FUNCTION hmpprt_grouplet_get_codelet_n(grouplet,codelet,name) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_GROUPLET_T), INTENT(IN)  :: grouplet  !< ...
       TYPE(HMPPRT_CODELET_T) , INTENT(OUT) :: codelet  !< ...
       CHARACTER(*)              , INTENT(IN)  :: name  !< ...
     END FUNCTION hmpprt_grouplet_get_codelet_n
  END_DUMMY_INTERFACE

  ! ================================================================
  ! ==  HMPPRT_DEVICE_T
  ! ================================================================

  DUMMY_INTERFACE
     !>
     !> Get the first n devices of the specified target
     !>
     !> The returned value is the number of devices stored in
     !> the argument devices (and so not greater than n)
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_get_devices(target,n,devices) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                 :: out
       INTEGER(HMPPRT_TARGET_T), INTENT(IN)  :: target  !< ...
       INTEGER(HMPPRT_INT_T),    INTENT(IN)  :: n        !< ...
       TYPE(HMPPRT_DEVICE_T), INTENT(OUT) :: devices(n)  !< ...
     END FUNCTION hmpprt_get_devices
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> \brief Return the number of available devices of the given target
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_get_device_count(target) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                :: out
       INTEGER(HMPPRT_TARGET_T), INTENT(IN) :: target  !< is any target including REF(hmpprt_enum,HMPPRT_HOST_TARGET)
     END FUNCTION hmpprt_get_device_count
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Get the device of a given target at a given rank
     !>
     !> \return HMPPRT_SUCCESS in case of success
     !>
     !> \remark Numbering starts at 1 in the Fortran binding
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_get_device(device,target,rank) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)               :: out
       TYPE(HMPPRT_DEVICE_T), INTENT(OUT)    :: device !< is the result
       INTEGER(HMPPRT_TARGET_T), INTENT(IN)  :: target !< can be any target including HMPPRT_ANY_TARGET
       INTEGER(HMPPRT_INT_T),    INTENT(IN)  :: rank   !< is the rank of the requested device (within the specified target)
     END FUNCTION hmpprt_get_device
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> Provide a readable name describing the device (e.g. "GeForce GTS 450")
     !>
     !> \ingroup grp_device
     SUBROUTINE hmpprt_device_get_name(device,name)
       USE hmpprt_enum
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device  !< is any device
       CHARACTER(*)          , INTENT(OUT) :: name    !< is the result
     END SUBROUTINE hmpprt_device_get_name
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Acquire a device.
     !>
     !> The following return codes can be expected:
     !>   - \c HMPPRT_SUCCESS        in case of success
     !>   - \c HMPPRT_LOCKFAIL_ERROR if this device is already in use.
     !>   - \c HMPPRT_IO_ERROR       if an error occurred while communicating with the resource manager.
     !>   - \c HMPPRT_UNKNOWN_ERROR  if an unexpected error occured.
     !>
     !> \see hmpprt_device_try_acquire
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_device_acquire(device) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)               :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN) :: device  !< ...
     END FUNCTION hmpprt_device_acquire
  END_DUMMY_INTERFACE


#ifdef DOC_MODE
  !> \brief Attempt to acquire a device.
  !>
  !> \return .TRUE. in case of success else .FALSE.
  !>
  !> \see hmpprt_device_acquire
  !>
  !> \ingroup grp_device
  FUNCTION hmpprt_device_try_acquire(device) result(out)
    USE hmpprt_enum
    BOOL_TYPE                          :: out
    TYPE(HMPPRT_DEVICE_T) , INTENT(IN) :: device  !< ...
  END FUNCTION hmpprt_device_try_acquire
#else

  INTERFACE hmpprt_device_try_acquire

#ifdef TRUE_IS_PLUS_1
     FUNCTION hmpprt_device_try_acquire_p1(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                          :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN) :: device  !< ...
     END FUNCTION hmpprt_device_try_acquire_p1
#endif


#ifdef TRUE_IS_MINUS_1
     FUNCTION hmpprt_device_try_acquire_m1(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                          :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN) :: device  !< ...
     END FUNCTION hmpprt_device_try_acquire_m1
#endif

  END INTERFACE hmpprt_device_try_acquire
#endif


  DUMMY_INTERFACE
     !>
     !> Release a previously acquired device.
     !>
     !> \ingroup grp_device
    FUNCTION hmpprt_device_release(device) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)               :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)    :: device  !< is any previously acquired device
     END FUNCTION hmpprt_device_release
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Perform a synchronous call of the specified  codelet with
     !> the specified argument list  args on the specified  device
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_device_call(device,codelet,args) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)             :: out
       TYPE(HMPPRT_DEVICE_T)  , INTENT(IN) :: device   !< is any device
       TYPE(HMPPRT_CODELET_T) , INTENT(IN) :: codelet  !< is any codelet compatible with that device
       TYPE(HMPPRT_ARGLIST_T) , INTENT(IN) :: args     !< is the codelet argument list
     END FUNCTION hmpprt_device_call
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> \brief Provide the default memory space for this device (usually the
     !> so-called global or main memory of the device)
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_device_get_default_memspace(device) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_MEMSPACE_T)             :: out
       TYPE(HMPPRT_DEVICE_T)  , INTENT(IN)    :: device  !< is any device
     END FUNCTION hmpprt_device_get_default_memspace
  END_DUMMY_INTERFACE


#ifdef DOC_MODE
  !>
  !> hmpprt_device_match_target(device,target)
  !>
  !> Return .TRUE. if the device match the specified target else .FALSE.
  !>
  !> \ingroup grp_device
  FUNCTION hmpprt_device_match_target(device,target) result(out)
    USE hmpprt_enum
    BOOL_TYPE                             :: out
    TYPE(HMPPRT_DEVICE_T)    , INTENT(IN) :: device   !< is any device
    INTEGER(HMPPRT_TARGET_T) , INTENT(IN) :: target   !< is the target to test
  END FUNCTION hmpprt_device_match_target
#else

  INTERFACE hmpprt_device_match_target

#ifdef TRUE_IS_PLUS_1
     FUNCTION hmpprt_device_match_target_p1(device,target) result(out)
       USE hmpprt_enum
       BOOL_TYPE                             :: out
       TYPE(HMPPRT_DEVICE_T) ,    INTENT(IN) :: device
       INTEGER(HMPPRT_TARGET_T) , INTENT(IN) :: target
     END FUNCTION hmpprt_device_match_target_p1
#endif

#ifdef TRUE_IS_MINUS_1
     FUNCTION hmpprt_device_match_target_m1(device,target) result(out)
       USE hmpprt_enum
       BOOL_TYPE                             :: out
       TYPE(HMPPRT_DEVICE_T) ,    INTENT(IN) :: device
       INTEGER(HMPPRT_TARGET_T) , INTENT(IN) :: target
     END FUNCTION hmpprt_device_match_target_m1
#endif
  END INTERFACE hmpprt_device_match_target
#endif



  DUMMY_INTERFACE
     !>
     !> Provide a readable name describing the device vendor (e.g. "NVIDIA")
     !>
     !> \ingroup grp_device
     SUBROUTINE hmpprt_device_get_vendor(device,vendor)
       USE hmpprt_enum
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device  !< is any device
       CHARACTER(*)          , INTENT(OUT) :: vendor  !< is the result
     END SUBROUTINE hmpprt_device_get_vendor
  END_DUMMY_INTERFACE


  !>
  !> Return .TRUE. if the device supports DOUBLE PRECISION
  !> else return .FALSE.
  !>
  !> \ingroup grp_device
#ifdef DOC_MODE
  FUNCTION hmpprt_device_double_support(device) result(out)
    USE hmpprt_enum
    BOOL_TYPE                           :: out
    TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device  !< is any device
  END FUNCTION hmpprt_device_double_support
#else

  INTERFACE hmpprt_device_double_support

#ifdef TRUE_IS_PLUS_1
     FUNCTION hmpprt_device_double_support_p1(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                           :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device
     END FUNCTION hmpprt_device_double_support_p1
#endif

#ifdef TRUE_IS_MINUS_1
     FUNCTION hmpprt_device_double_support_m1(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                           :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device
     END FUNCTION hmpprt_device_double_support_m1
#endif

  END INTERFACE hmpprt_device_double_support

#endif


  !>
  !> Return .TRUE. if direct access to the host memory is supported by this device
  !> else return .FALSE.
  !>
  !> \ingroup grp_device
#ifdef DOC_MODE
  FUNCTION hmpprt_device_hostmem_support(device) result(out)
    USE hmpprt_enum
    BOOL_TYPE                           :: out
    TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device  !< is any device
  END FUNCTION hmpprt_device_hostmem_support
#else

  INTERFACE hmpprt_device_hostmem_support

#ifdef TRUE_IS_PLUS_1
     FUNCTION hmpprt_device_hostmem_support_p1(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                           :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device
     END FUNCTION hmpprt_device_hostmem_support_p1
#endif

#ifdef TRUE_IS_MINUS_1
     FUNCTION hmpprt_device_hostmem_support_m1(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                           :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device
     END FUNCTION hmpprt_device_hostmem_support_m1
#endif

  END INTERFACE hmpprt_device_hostmem_support

#endif


  DUMMY_INTERFACE
     !>
     !> Provide the max clock rate of the device expressed in Mhz
     !>
     !> \return 0 if the device does not provide this information
     !>
     !> \ingroup grp_device
    FUNCTION hmpprt_device_get_max_clock_rate(device) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)               :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device      !< is any device
     END FUNCTION hmpprt_device_get_max_clock_rate
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Allocate size bytes of memory in the given memory space of a device.
     !>
     !> The result is device and memory space dependant but it can be safely
     !> assumed that a return value of 0 indicates a failure.
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_device_allocate(data) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)          :: out
       TYPE(HMPPRT_DATA_T) , INTENT(IN) :: data !< ...
     END FUNCTION hmpprt_device_allocate
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> Free a memory area previously allocated by hmpprt_device_allocate
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_device_free(data) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)          :: out
       TYPE(HMPPRT_DATA_T) , INTENT(IN) :: data !< ...
     END FUNCTION hmpprt_device_free
  END_DUMMY_INTERFACE

#ifdef DOC_MODE
  DUMMY_INTERFACE
     !>
     !> Upload size bytes from the specified host address to a device address
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_device_upload(data,from,offset,size) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_DATA_T)        , INTENT(IN) :: data     !< is the destination data
       CLASS(*)                   , INTENT(IN) :: from(*)  !< is the source (on the host)
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN) :: offset   !< is the offset in bytes
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN) :: size     !< is the amount of bytes to upload
     END FUNCTION hmpprt_device_upload
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE
     FUNCTION hmpprt_device_upload(data,from,offset,size) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_DATA_T)        , INTENT(IN) :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK :: from
       INTEGER                    , INTENT(IN) :: from(*)
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN) :: offset
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN) :: size
     END FUNCTION hmpprt_device_upload
  END_DUMMY_INTERFACE
#else
  INTEGER(HMPPRT_ERROR_T) , EXTERNAL :: hmpprt_device_upload
#endif

#ifdef DOC_MODE
  DUMMY_INTERFACE
     !>
     !> Download bytes from device to host
     !>
     !> \ingroup grp_device
     FUNCTION hmpprt_device_download(data,to,offset,size) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                    :: out
       TYPE(HMPPRT_DATA_T)        , INTENT(IN)    :: data     !< is the source data
       CLASS(*)                   , INTENT(INOUT) :: to(*)    !< is the destination (on the host)
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN)    :: offset   !< is the offset in bytes
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN)    :: size     !< is the amount of bytes to download

     END FUNCTION hmpprt_device_download
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE

     FUNCTION hmpprt_device_download(data,to,offset,size) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                    :: out
       TYPE(HMPPRT_DATA_T)        , INTENT(IN)    :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK :: to
       INTEGER                    , INTENT(INOUT) :: to(*)
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN)    :: offset
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN)    :: size
     END FUNCTION hmpprt_device_download
  END_DUMMY_INTERFACE
#else
  INTEGER(HMPPRT_ERROR_T) , EXTERNAL :: hmpprt_device_download
#endif


  ! ================================================================
  ! ==  HOST Device Support
  ! ================================================================

  DUMMY_INTERFACE
     !>
     !> Provide the first and unique HOST device
     !>
     !> \remark This is a shortcut for \ref hmpprt_get_device "hmpprt_get_device(HMPPRT_HOST_DEVICE,1)"
     !>
     !> \ingroup grp_device
     SUBROUTINE hmpprt_get_host_device(device)
       USE hmpprt_enum
       TYPE(HMPPRT_DEVICE_T), INTENT(OUT) :: device !< is the result
     END SUBROUTINE  hmpprt_get_host_device
  END_DUMMY_INTERFACE

  ! ================================================================
  ! ==  CUDA Device Support
  ! ================================================================

  DUMMY_INTERFACE
     !>
     !> Get the first \c n CUDA devices
     !>
     !> \return The number of devices effectively stored in \c devices. That
     !> value cannot be greater than \c n.
     !>
     !> \ingroup grp_cuda
     FUNCTION hmpprt_cuda_get_devices(n,devices) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                 :: out
       INTEGER(HMPPRT_INT_T),  INTENT(IN)    :: n           !< is the size of the argument \c devices
       TYPE(HMPPRT_DEVICE_T),  INTENT(OUT)   :: devices(n)  !< is the result
     END FUNCTION hmpprt_cuda_get_devices
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Return the number of available CUDA devices
     !>
     !> \remark This is a shortcut for hmpprt_get_device_count(HMPPRT_CUDA_TARGET).
     !>
     !> \ingroup grp_cuda
     FUNCTION hmpprt_cuda_get_device_count() result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                :: out   !< ...
     END FUNCTION hmpprt_cuda_get_device_count
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Get the CUDA device at the specified rank
     !>
     !> \return Return HMPPRT_SUCCESS in case of success
     !>
     !> \remark Numbering of devices starts at 1 in the Fortran API
     !>
     !> \ingroup grp_cuda
     FUNCTION hmpprt_cuda_get_device(device,rank) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)               :: out
       TYPE(HMPPRT_DEVICE_T),    INTENT(OUT) :: device   !< should be a CUDA device
       INTEGER(HMPPRT_INT_T),    INTENT(IN)  :: rank     !< is the rank of the requested device
     END FUNCTION hmpprt_cuda_get_device
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> \brief Return the version of the CUDA driver used by a CUDA device
     !>
     !> \ingroup grp_cuda
    FUNCTION hmpprt_cuda_get_driver_version(device) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)               :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device !< should be a CUDA device
     END FUNCTION hmpprt_cuda_get_driver_version
  END_DUMMY_INTERFACE



  DUMMY_INTERFACE
     !>
     !> \brief Provide the Compute Capability of a CUDA device
     !>
     !> \remark In case of failure (e.g. not a CUDA device), the arguments \c major
     !> and \c minor are both set to 0
     !>
     !> \ingroup grp_cuda
     SUBROUTINE hmpprt_cuda_get_compute_capability(device,major,minor)
       USE hmpprt_enum
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device  !< should be a CUDA device
       INTEGER(HMPPRT_INT_T) , INTENT(OUT) :: major   !< is set to the major part of the device compute capability
       INTEGER(HMPPRT_INT_T) , INTENT(OUT) :: minor   !< is set to the minor part of the device compute capability
     END SUBROUTINE  hmpprt_cuda_get_compute_capability
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> \brief Provide the amount of main memory on a CUDA device
     !>
     !> \return a valid size or 0 in case of failure
     !>
     !> \ingroup grp_cuda
     FUNCTION hmpprt_cuda_get_memory_size(device) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_SIZE_T)              :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device  !< should be a CUDA device
     END FUNCTION hmpprt_cuda_get_memory_size
  END_DUMMY_INTERFACE


  !> \return .TRUE. if memory ECC (error correction code) is enabled
  !> for a CUDA device else .FALSE.
  !>
  !> \ingroup grp_cuda
#ifdef DOC_MODE
   FUNCTION hmpprt_cuda_has_ecc(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                           :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device !< should be a CUDA device
     END FUNCTION hmpprt_cuda_has_ecc
#else
  INTERFACE hmpprt_cuda_has_ecc

#ifdef TRUE_IS_PLUS_1
     FUNCTION hmpprt_cuda_has_ecc_p1(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                           :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device
     END FUNCTION hmpprt_cuda_has_ecc_p1
#endif

#ifdef TRUE_IS_MINUS_1
     FUNCTION hmpprt_cuda_has_ecc_m1(device) result(out)
       USE hmpprt_enum
       BOOL_TYPE                           :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device
     END FUNCTION hmpprt_cuda_has_ecc_m1
#endif

  END INTERFACE hmpprt_cuda_has_ecc
#endif

  ! ================================================================
  ! ==  OPENCL Device Support
  ! ================================================================



  DUMMY_INTERFACE
     !>
     !> Get the first n OpenCL devices
     !>
     !> \return The number of devices effectively stored in \c devices. That
     !> value cannot be greater than \c n.
     !>
     !> \ingroup grp_opencl
     FUNCTION hmpprt_opencl_get_devices(n,devices) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                 :: out
       INTEGER(HMPPRT_INT_T),    INTENT(IN)  :: n          !< is the size of the argument \c devices
       TYPE(HMPPRT_DEVICE_T), INTENT(OUT)    :: devices(n) !< is the result
     END FUNCTION hmpprt_opencl_get_devices
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Return the number of available OpenCL devices
     !>
     !> \remark This is a shortcut for hmpprt_get_device_count(HMPPRT_OPENCL_TARGET).
     !>
     !> \ingroup grp_opencl
     FUNCTION hmpprt_opencl_get_device_count() result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                :: out
     END FUNCTION hmpprt_opencl_get_device_count
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Get the OpenCL device at the specified rank
     !>
     !> \return HMPPRT_SUCCESS in case of success
     !>
     !> \remark Numbering of devices starts at 1 in the Fortran API
     !>
     !> \ingroup grp_opencl
     FUNCTION hmpprt_opencl_get_device(device,rank) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)               :: out
       TYPE(HMPPRT_DEVICE_T),   INTENT(OUT)  :: device  !< is the result
       INTEGER(HMPPRT_INT_T),    INTENT(IN)  :: rank    !< is the rank of the requested device
     END FUNCTION hmpprt_opencl_get_device
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Provide the OpenCL driver version as a string
     !>
     !> \ingroup grp_opencl
     SUBROUTINE hmpprt_opencl_get_driver_version(device,version)
       USE hmpprt_enum
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)  :: device  !< should be an OpenCL device
       CHARACTER(*)          , INTENT(OUT) :: version !< is the result
     END SUBROUTINE  hmpprt_opencl_get_driver_version

  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Provide OpenCL device version as a string
     !>
     !> \ingroup grp_opencl
     SUBROUTINE hmpprt_opencl_get_device_version(device,version)
       USE hmpprt_enum
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)   :: device  !< should be an OpenCL device
       CHARACTER(*)          , INTENT(OUT) :: version  !< is the result
     END SUBROUTINE  hmpprt_opencl_get_device_version

  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> brief Provide the amount of main memory on a OpenCL device
     !>
     !> \return a valid size or 0 in case of failure
     !>
     !> \ingroup grp_opencl
     FUNCTION hmpprt_opencl_get_memory_size(device) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_SIZE_T)                 :: out
       TYPE(HMPPRT_DEVICE_T) , INTENT(IN)     :: device  !< should be an OpenCL device
     END FUNCTION hmpprt_opencl_get_memory_size
  END_DUMMY_INTERFACE



  ! ================================================================
  ! ==  HMPPRT_DATA_T
  ! ================================================================


  DUMMY_INTERFACE
     !> Create a device data descriptor of the specified size (in bytes) in the
     !> default memory space of the device.
     !>
     !> The data itself is not yet allocated (see \ref hmpprt_data_allocate_content)
     !>
     !> This function is shortcut for
     !>   hmpprt_data_new_m(device, size, hmpprt_device_get_default_memspace(device))
     !>
      !> \ingroup grp_data
    SUBROUTINE hmpprt_data_new(data, device, size)
       USE hmpprt_enum
       TYPE(HMPPRT_DATA_T)     , INTENT(OUT)  :: data  !< ...
       TYPE(HMPPRT_DEVICE_T)   , INTENT(IN) :: device  !< ...
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN) :: size !< ...
     END SUBROUTINE  hmpprt_data_new
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !> Create a device data descriptor of the specified size (in bytes) in the
     !> specified memory space of the device.
     !>
     !> The data itself is not yet allocated (see \ref hmpprt_data_allocate_content)
     !>
     !> \ingroup grp_data
     SUBROUTINE hmpprt_data_new_m(data, device, size, memspace)
       USE hmpprt_enum

       TYPE(HMPPRT_DATA_T)        , INTENT(OUT) :: data      !< ...
       TYPE(HMPPRT_DEVICE_T)      , INTENT(IN)  :: device    !< The device on which the data will reside
       INTEGER(HMPPRT_SIZE_T)     , INTENT(IN)  :: size      !< The data size in bytes
       INTEGER(HMPPRT_MEMSPACE_T) , INTENT(IN)  :: memspace  !< The memory space on the device
     END SUBROUTINE hmpprt_data_new_m

     !>
     !> Delete a data descriptor previously created with
     !> \ref hmpprt_data_new or \ref hmpprt_data_new_m
     !>
     !> \ingroup grp_data
     SUBROUTINE hmpprt_data_delete(data)
       USE hmpprt_enum
       TYPE(HMPPRT_DATA_T)     , INTENT(IN) :: data     !< ...
     END SUBROUTINE hmpprt_data_delete
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> Provide the device on which a data descriptor operates
     !>
      !> \ingroup grp_data
    SUBROUTINE hmpprt_data_get_device(data,device)
       USE hmpprt_enum
       TYPE(HMPPRT_DATA_T)     , INTENT(IN)  :: data    !< ...
       TYPE(HMPPRT_DEVICE_T)   , INTENT(OUT) :: device  !< ...
     END SUBROUTINE  hmpprt_data_get_device
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> Provide the memory space on which a data descriptor operates
     !>
     !> \ingroup grp_data
     FUNCTION hmpprt_data_get_memspace(data) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_MEMSPACE_T)           :: out
       TYPE(HMPPRT_DATA_T)     , INTENT(IN) :: data !< ...
     END FUNCTION hmpprt_data_get_memspace
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Provide the device address currently managed by this data descriptor
     !>
     !> \ingroup grp_data
     SUBROUTINE hmpprt_data_get_device_address(data,address)
       USE hmpprt_enum
       INTEGER(HMPPRT_ADDRESS_T)  , INTENT(OUT) :: address   !< ...
       TYPE(HMPPRT_DATA_T)        , INTENT(IN)  :: data   !< ...
     END SUBROUTINE hmpprt_data_get_device_address
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Set the device address to be managed by this data descriptor
     !>
     !> Remark: This function is for advanced users only.
     !>         The recommanded way to associate a device address to a
     !>         data descriptor is via an allocation function such
     !>         as hmpprt_data_allocate_content()
     !>
     !> \ingroup grp_data
     SUBROUTINE hmpprt_data_set_device_address(data,address)
       USE hmpprt_enum
       TYPE(HMPPRT_DATA_T)     , INTENT(IN)  :: data   !< ...
       INTEGER(HMPPRT_ADDRESS_T)  , INTENT(IN) :: address   !< ...
     END SUBROUTINE hmpprt_data_set_device_address
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !> Allocate memory in a data descriptor using the size
     !> specified at its creation
      !> \ingroup grp_data
    FUNCTION hmpprt_data_allocate_content(data) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_DATA_T)     , INTENT(IN) :: data   !< ...
     END FUNCTION hmpprt_data_allocate_content
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !> Free the memory previously allocated using a function
     !> such as hmpprt_data_allocate_content()
     !> \ingroup grp_data
     FUNCTION hmpprt_data_free_content(data) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                 :: out
       TYPE(HMPPRT_DATA_T)     , INTENT(IN) :: data   !< ...
     END FUNCTION hmpprt_data_free_content
  END_DUMMY_INTERFACE

#ifdef DOC_MODE
  DUMMY_INTERFACE
     !> Upload the whole data from host to device memory
     !>
     !> \ingroup grp_data
     FUNCTION hmpprt_data_upload(data,host_address) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                    :: out
       TYPE(HMPPRT_DATA_T),           INTENT(IN)  :: data            !< ...
       CLASS(*),                      INTENT(IN)  :: host_address(*) !< ...
     END FUNCTION hmpprt_data_upload
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
   DUMMY_INTERFACE
     !> Upload the whole data from host to device memory
     !>
     !> \ingroup grp_data
     FUNCTION hmpprt_data_upload(data,host_address) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                    :: out
       TYPE(HMPPRT_DATA_T),        INTENT(IN)     :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  host_address
       INTEGER,                    INTENT(IN)     :: host_address(*)
     END FUNCTION hmpprt_data_upload
  END_DUMMY_INTERFACE
#else
  INTEGER(HMPPRT_ERROR_T) , EXTERNAL :: hmpprt_data_upload
#endif

#ifdef DOC_MODE
  DUMMY_INTERFACE
     !> Upload part of the data from host to device memory
     !>
     !> \ingroup grp_data
     FUNCTION hmpprt_data_upload_part(data,host_address,offset,size) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                    :: out
       TYPE(HMPPRT_DATA_T),           INTENT(IN)  :: data
       CLASS(*),                      INTENT(IN)  :: host_address(*)
       INTEGER(HMPPRT_SIZE_T),        INTENT(IN)  :: offset
       INTEGER(HMPPRT_SIZE_T),        INTENT(IN)  :: size
     END FUNCTION hmpprt_data_upload_part
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE
     !> Upload part of the data from host to device memory
     !>
     !> \ingroup grp_data
     FUNCTION hmpprt_data_upload_part(data,host_address,offset,size) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                    :: out
       TYPE(HMPPRT_DATA_T),           INTENT(IN)  :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  host_address
       INTEGER,                       INTENT(IN)  :: host_address(*)
       INTEGER(HMPPRT_SIZE_T),        INTENT(IN)  :: offset
       INTEGER(HMPPRT_SIZE_T),        INTENT(IN)  :: size
     END FUNCTION hmpprt_data_upload_part
  END_DUMMY_INTERFACE
#else
  INTEGER(HMPPRT_ERROR_T) , EXTERNAL :: hmpprt_data_upload_part
#endif

#ifdef DOC_MODE
  DUMMY_INTERFACE
     !> Download the whole data from device to host memory
     !>
     !> \ingroup grp_data
     FUNCTION hmpprt_data_download(data,host_address) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                      :: out
       TYPE(HMPPRT_DATA_T),        INTENT(IN)       :: data              !< ...
       CLASS(*),                   INTENT(INOUT)    :: host_address(*)   !< ...
     END FUNCTION hmpprt_data_download
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE
     FUNCTION hmpprt_data_download(data,host_address) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                      :: out
       TYPE(HMPPRT_DATA_T),        INTENT(IN)       :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  host_address
       INTEGER,                    INTENT(INOUT)    :: host_address(*)
     END FUNCTION hmpprt_data_download
  END_DUMMY_INTERFACE
#else
  INTEGER(HMPPRT_ERROR_T) , EXTERNAL :: hmpprt_data_download
#endif


#ifdef DOC_MODE
  DUMMY_INTERFACE
     !> Download the part of the data from device to host memory
     !>
     !> \ingroup grp_data
     FUNCTION hmpprt_data_download_part(data,host_address,offset,size) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                   :: out
       TYPE(HMPPRT_DATA_T),        INTENT(IN)    :: data              !< ...
       CLASS(*),                   INTENT(INOUT) :: host_address(*)   !< ...
       INTEGER(HMPPRT_SIZE_T),     INTENT(IN)    :: offset
       INTEGER(HMPPRT_SIZE_T),     INTENT(IN)    :: size
     END FUNCTION hmpprt_data_download_part
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE
     FUNCTION hmpprt_data_download_part(data,host_address,offset,size) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                   :: out
       TYPE(HMPPRT_DATA_T),        INTENT(IN)    :: data              !< ...
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  host_address
       INTEGER,                    INTENT(INOUT) :: host_address(*)   !< ...
       INTEGER(HMPPRT_SIZE_T),     INTENT(IN)    :: offset
       INTEGER(HMPPRT_SIZE_T),     INTENT(IN)    :: size
     END FUNCTION hmpprt_data_download_part
  END_DUMMY_INTERFACE
#else
  INTEGER(HMPPRT_ERROR_T) , EXTERNAL :: hmpprt_data_download_part
#endif


  ! ================================================================
  ! ==  HMPPRT_ARGLIST_T
  ! ================================================================


  DUMMY_INTERFACE
     !>
     !> Create an empty argument list
     !>
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_new(args)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T)  , INTENT(OUT) :: args   !< is the new argument list
     END SUBROUTINE hmpprt_arglist_new
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Create an exact copy of an argument list
     !>
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_clone(args_new,args_old)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(OUT) :: args_new   !< is the new argument list
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args_old   !< is the old argument list
     END SUBROUTINE hmpprt_arglist_clone
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Provide the number of elements currently stored in the argument list
     !>
     !> \ingroup grp_arglist
     FUNCTION hmpprt_arglist_size(args) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                  :: out
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is an argument list
     END FUNCTION hmpprt_arglist_size
  END_DUMMY_INTERFACE

  !>
  !>
  !>
  DUMMY_INTERFACE
     !>
     !> \brief Add an argument whose value is provided by a Data buffer descriptor
     !>
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_data(args,data)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       TYPE(HMPPRT_DATA_T),       INTENT(IN)  :: data   !< is the Data buffer
     END SUBROUTINE hmpprt_arglist_add_data

  END_DUMMY_INTERFACE

#ifdef DOC_MODE
  DUMMY_INTERFACE
     !>
     !> Add arbitrary data specified by a start address and a number of bytes
     !>
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_value(args,value,size)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args     !< is the argument list to be extended
       CLASS(*) ,                 INTENT(IN)  :: value(*) !< is any object whose value must be passed
       INTEGER(HMPPRT_SIZE_T) ,   INTENT(IN)  :: size     !< is the size of data in bytes
     END SUBROUTINE hmpprt_arglist_add_value
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE
     SUBROUTINE hmpprt_arglist_add_value(args,value,size)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  value
       INTEGER,                   INTENT(IN)  :: value(*)
       INTEGER(HMPPRT_SIZE_T) ,   INTENT(IN)  :: size
     END SUBROUTINE hmpprt_arglist_add_value
  END_DUMMY_INTERFACE
#else
    EXTERNAL hmpprt_arglist_add_value
#endif



  !
  ! REMINDER: If a function is added to this interface then a
  !           MODULE PROCEDURE entry should be added to the interface
  !           above (for Doxygen)
  !
  NAMED_INTERFACE(hmpprt_arglist_add_scalar)

#ifdef KIND_INT8

     !> Add an argument of type 8bit integer to an argument list
     !> \see The generic interface hmpprt_arglist_add_scalar
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_int8(args,value)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       INTEGER(KIND_INT8),        INTENT(IN)  :: value  !< is the value to be added
     END SUBROUTINE hmpprt_arglist_add_int8

#endif

#ifdef KIND_INT16

     !> Add an argument of type 16bit integer to an argument list
     !> \see The generic interface hmpprt_arglist_add_scalar
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_int16(args,value)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       INTEGER(KIND_INT16),       INTENT(IN)  :: value  !< is the value to be added
     END SUBROUTINE hmpprt_arglist_add_int16

#endif

#ifdef KIND_INT32

     !> Add an argument of type 32bit integer to an argument list
     !> \see The generic interface hmpprt_arglist_add_scalar
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_int32(args,value)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       INTEGER(KIND_INT32),       INTENT(IN)  :: value  !< is the value to be added
     END SUBROUTINE hmpprt_arglist_add_int32

#endif

#ifdef KIND_INT64

     !> Add an argument of type 64bit integer to an argument list
     !> \see The generic interface hmpprt_arglist_add_scalar
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_int64(args,value)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       INTEGER(KIND_INT64),       INTENT(IN)  :: value  !< is the value to be added
     END SUBROUTINE hmpprt_arglist_add_int64

#endif

#ifdef KIND_REAL32

     !> Add an argument of type 32bit floating point to an argument list
     !> \see The generic interface hmpprt_arglist_add_scalar
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_real32(args,value)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       REAL(KIND_REAL32),         INTENT(IN)  :: value  !< is the value to be added
     END SUBROUTINE hmpprt_arglist_add_real32

     !> Add an argument of type 2x32bit complex to an argument list
     !> \see The generic interface hmpprt_arglist_add_scalar
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_complex32(args,value)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       COMPLEX(KIND_REAL32),      INTENT(IN)  :: value  !< is the value to be added
     END SUBROUTINE hmpprt_arglist_add_complex32

#endif

#ifdef KIND_REAL64

     !> Add an argument of type 64bit floating point to an argument list
     !> \see The generic interface hmpprt_arglist_add_scalar
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_real64(args,value)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       REAL(KIND_REAL64),         INTENT(IN)  :: value  !< is the value to be added
     END SUBROUTINE hmpprt_arglist_add_real64

     !> Add an argument of type 2x64bit complex to an argument list
     !> \see The generic interface hmpprt_arglist_add_scalar
     !> \ingroup grp_arglist
     SUBROUTINE hmpprt_arglist_add_complex64(args,value)
       USE hmpprt_enum
       TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args   !< is the argument list to be extended
       COMPLEX(KIND_REAL64),      INTENT(IN)  :: value  !< is the value to be added
     END SUBROUTINE hmpprt_arglist_add_complex64

#endif

  END_NAMED_INTERFACE(hmpprt_arglist_add_scalar)

  !> Add multiple arguments to an arglist using a format.
  !>
  !> The number and type of each argument is specified by a format string
  !> in which each character (after the leading '(' and the trailing ')')
  !> represents an argument:
  !>
  !>    - \c 'B' or \c 'b' for a 1 byte  integer scalar (or any other type of that size, usually \c INTEGER(1))
  !>    - \c 'W' or \c 'w' for a 2 bytes integer scalar (or any other type of that size, usually \c INTEGER(2))
  !>    - \c 'I' or \c 'i' for a 4 bytes integer scalar (or any other type of that size, usually \c INTEGER(4))
  !>    - \c 'L' or \c 'l' for a 8 bytes integer scalar (or any other type of that size, usually \c INTEGER(8))
  !>    - \c 'R' or \c 'r' for a 32 bit \c REAL    scalar (single precision, usually \c REAL(4))
  !>    - \c 'D' or \c 'd' for a 64 bit \c REAL    scalar (double precision, usually \c REAL(8))
  !>    - \c 'C' or \c 'c' for a 32 bit \c COMPLEX scalar (single precision, usually \c COMPLEX(4))
  !>    - \c 'Z' or \c 'z' for a 64 bit \c COMPLEX scalar (double precision, usually \c COMPLEX(8) )
  !>    - \c '@'           for a data descriptor of type \c TYPE(REF(hmpprt_enum,HMPPRT_DATA_T))
  !>    - \c 'M' or \c 'm' for a memory area specified by an arbitrary fortran expression and an \c INTEGER(REF(hmpprt_enum,HMPPRT_SIZE_T) )
  !>                       giving a size in bytes.
  !>    - \c '#'           to insert a marker of value REF(hmpprt_enum,HMPPRT_ARG_CHECK). The value of that marker
  !>                       is checked at runtime by HMPP
  !>
  !>
  !> Here is an example showing how an argument list of length 7
  !> can be easily constructed
  !>
  !> \verbatim
  !>     TYPE(HMPPRT_ARGLIST_T)    :: args
  !>     TYPE(HMPPRT_DATA_T)       :: data_A , data_B
  !>     INTEGER                   :: na, nb      !> some 32bit integers
  !>     REAL(4)                   :: pi          !> a 32 bit real
  !>     REAL(8)                   :: omega       !> a 64bit real
  !>     CHARACTER(20)             :: name
  !>     ...
  !>     call hmpprt_arglist_fadd(args, "(ii@@rdm#)" &
  !>       & , na, nb                 &
  !>       & , data_A, data_B         &
  !>       & , pi, omega &            &
  !>       & , name, 20_HMPPRT_SIZE_T &
  !>       & , HMPPRT_ARG_CHECK       &
  !>       & )
  !> \endverbatim
  !>
  !> \ingroup grp_arglist
#ifdef VARARG_SUPPORT
  SUBROUTINE hmpprt_arglist_fadd(args,format,__vararg__)
    USE hmpprt_enum
    TYPE(HMPPRT_ARGLIST_T),    INTENT(IN)  :: args        !< is the argument list to be extended
    CHARACTER(*),              INTENT(IN)  :: format      !< is a string format (see above)
    CLASS(*)                   INTENT(IN)  :: __vararg__  !< any number of arguments
  END SUBROUTINE hmpprt_arglist_fadd
#else
  EXTERNAL hmpprt_arglist_fadd
#endif
  ! ================================================================
  ! ==  HMPPRT_CODELET_T
  ! ================================================================

  DUMMY_INTERFACE
     !> \brief Provide the grouplet that owns a codelet
     !>
     !> \ingroup grp_codelet
     SUBROUTINE hmpprt_codelet_get_grouplet(codelet,grouplet)
       USE hmpprt_enum
       TYPE(HMPPRT_CODELET_T),  INTENT(IN)   :: codelet   !< is the considered codelet
       TYPE(HMPPRT_GROUPLET_T), INTENT(OUT)  :: grouplet  !< is the result
     END SUBROUTINE hmpprt_codelet_get_grouplet
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Provide the memory space of an argument specified by its rank
     !>
     !> \remark In the HMPPRT Fortran binding, arguments numbering starts at 1
     !> \ingroup grp_codelet
     !>
     FUNCTION hmpprt_codelet_get_arg_memspace_i(codelet,rank) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_MEMSPACE_T)              :: out
       TYPE(HMPPRT_CODELET_T),    INTENT(IN)   :: codelet  !< is the considered codelet
       INTEGER(HMPPRT_INT_T)    , INTENT(IN)   :: rank     !< is the rank of the argument
     END FUNCTION hmpprt_codelet_get_arg_memspace_i
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !> \brief Provide the memory space of an argument specified by its name.
     !>
     !> \ingroup grp_codelet
     FUNCTION hmpprt_codelet_get_arg_memspace_n(codelet,name) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_MEMSPACE_T)            :: out
       TYPE(HMPPRT_CODELET_T) , INTENT(IN)   :: codelet   !<
       CHARACTER(*)           , INTENT(IN)   :: name      !< is the name of the argument
     END FUNCTION hmpprt_codelet_get_arg_memspace_n
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !> \brief Provide the name of a codelet
     !>
     !> \ingroup grp_codelet
     SUBROUTINE hmpprt_codelet_get_name(codelet,name)
       USE hmpprt_enum
       TYPE(HMPPRT_CODELET_T), INTENT(IN)   :: codelet   !<
       CHARACTER(*)          , INTENT(OUT)  :: name      !< contains the codelet name on exit
     END SUBROUTINE hmpprt_codelet_get_name
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !> \brief Call a codelet on a given device.
     !>
     !> The codelet and its arguments must be compatible with that device
     !>
     !> \ingroup grp_codelet
     FUNCTION hmpprt_codelet_call(codelet, device, args) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)                :: out
       TYPE(HMPPRT_CODELET_T)   , INTENT(IN)  :: codelet   !< is the codelet to be call
       TYPE(HMPPRT_DEVICE_T)    , INTENT(IN)  :: device    !< is the device on which the call should occur
       TYPE(HMPPRT_ARGLIST_T)   , INTENT(IN)  :: args      !< is the argument list to use
     END FUNCTION hmpprt_codelet_call
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !> \brief Provide the number of arguments expected by a codelet
     !> \ingroup grp_codelet
     FUNCTION hmpprt_codelet_get_arg_count(codelet) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_INT_T)                :: out
       TYPE(HMPPRT_CODELET_T), INTENT(IN)   :: codelet   !< is any codelet
     END FUNCTION hmpprt_codelet_get_arg_count
  END_DUMMY_INTERFACE


  ! ================================================================
  ! ==  HMPPRT_QUEUE_T
  ! ================================================================

  !> \addtogroup grp_queues
  !> @{

  DUMMY_INTERFACE
     !>
     !> Create a new execution queue
     !>
     SUBROUTINE hmpprt_queue_new(queue)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(OUT)  :: queue   !< is any queue
     END SUBROUTINE  hmpprt_queue_new
  END_DUMMY_INTERFACE
  !> @}




  DUMMY_INTERFACE
     !>
     !> \brief  Delete an execution queue previously created with hmpprt_queue_new
     !> \ingroup grp_queues
     !>
     SUBROUTINE hmpprt_queue_delete(queue)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue   !< is any queue
     END SUBROUTINE hmpprt_queue_delete
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> \brief  Launch an asynchronous execution of the queue operations.
     !> \ingroup grp_queues
     !>
     !> Once the queue is running, it is still possible to enqueue an operation.
     !>
     SUBROUTINE hmpprt_queue_start(queue)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue   !< is any queue
     END SUBROUTINE hmpprt_queue_start
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> \brief Wait until every asynchronous operation in the queue is finished
     !> \ingroup grp_queues
     !>
     FUNCTION hmpprt_queue_wait(queue) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)             :: out
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue   !< is any queue
     END FUNCTION hmpprt_queue_wait
  END_DUMMY_INTERFACE

    !>   @}


  DUMMY_INTERFACE
     !>
     !> \brief Similar to hmpprt_queue_wait() but only wait for up to timeout_ms
     !> milliseconds.
     !>
     !> if the timeout is reached the result result is HMPPRT_TIMEOUT_ERROR
     !>
     !> HMPPRT_SUCCESS is returned if the queue was successfully executed
     !> else the error code produced by the failing operation is returned
     !>
     !> \ingroup grp_queues
     FUNCTION hmpprt_queue_wait_t(queue,timeout_ms) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)             :: out
       TYPE(HMPPRT_QUEUE_T),    INTENT(IN) :: queue      !< is any queue
       INTEGER(HMPPRT_SIZE_T),  INTENT(IN) :: timeout_ms !< is a timeout in milli-seconds
     END FUNCTION hmpprt_queue_wait_t
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> Perform a synchronous execution of the queue.
     !>
     !> \ingroup grp_queues
     FUNCTION hmpprt_queue_execute(queue) result(out)
       USE hmpprt_enum
       INTEGER(HMPPRT_ERROR_T)          :: out
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue   !< is any queue
     END FUNCTION hmpprt_queue_execute
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !> Clear the content of a queue, removing every operations without executing them.
     !>
     !> \ingroup grp_queues
     SUBROUTINE hmpprt_queue_clear(queue)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue   !< ...
     END SUBROUTINE hmpprt_queue_clear
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !>  Enqueue the execution of REFCALL(hmpprt_device_acquire,(device))
     !> \ingroup grp_queues
     !>
     SUBROUTINE hmpprt_enqueue_acquire(queue,device)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN)  :: queue  !< is any queue
       TYPE(HMPPRT_DEVICE_T), INTENT(IN) :: device !< is any device
     END SUBROUTINE hmpprt_enqueue_acquire
  END_DUMMY_INTERFACE



  DUMMY_INTERFACE
     !>
     !>  Enqueue the execution of REFCALL(hmpprt_device_call,(device,codelet,args))
     !>
     !> \ingroup grp_queues
     SUBROUTINE hmpprt_enqueue_call(queue,device,codelet,args)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T),  INTENT(IN)   :: queue    !< is any queue
       TYPE(HMPPRT_DEVICE_T), INTENT(IN)   :: device   !< is any device
       TYPE(HMPPRT_CODELET_T) , INTENT(IN) :: codelet  !< is any codelet
       TYPE(HMPPRT_ARGLIST_T) , INTENT(IN) :: args     !< is an argument list
     END SUBROUTINE hmpprt_enqueue_call
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !>  Enqueue the execution of REFCALL(hmpprt_data_allocate_content,(data))
     !>
     !> \ingroup grp_queues
     SUBROUTINE hmpprt_enqueue_allocate(queue,data)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue !< is any queue
       TYPE(HMPPRT_DATA_T) , INTENT(IN) :: data  !< is any data buffer descriptor
     END SUBROUTINE hmpprt_enqueue_allocate
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !>  Enqueue the execution of REFCALL(hmpprt_data_free_content,(data))
     !>
     !> \ingroup grp_queues
     SUBROUTINE hmpprt_enqueue_free(queue,data)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue !< is any queue
       TYPE(HMPPRT_DATA_T) , INTENT(IN) :: data  !< is any data buffer descriptor
     END SUBROUTINE hmpprt_enqueue_free
  END_DUMMY_INTERFACE




#ifdef DOC_MODE
  DUMMY_INTERFACE
     !>
     !> Enqueue the execution of REFCALL(hmpprt_data_upload,(data,host_address))
     !>
     !> \ingroup grp_queues
     SUBROUTINE hmpprt_enqueue_upload(queue,data,host_address)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN)    :: queue !< ...
       TYPE(HMPPRT_DATA_T),  INTENT(IN)    :: data !< ...
       CLASS(*),             INTENT(IN)    :: host_address(*) !< ...
     END SUBROUTINE hmpprt_enqueue_upload
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE
     SUBROUTINE hmpprt_enqueue_upload(queue,data,host_address)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN)    :: queue
       TYPE(HMPPRT_DATA_T),  INTENT(IN)    :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  host_address
       INTEGER,              INTENT(IN)    :: host_address(*)
     END SUBROUTINE hmpprt_enqueue_upload
  END_DUMMY_INTERFACE
#else
  EXTERNAL hmpprt_enqueue_upload
#endif


#ifdef DOC_MODE
  DUMMY_INTERFACE
     !>
     !> Enqueue the execution of REFCALL(hmpprt_data_upload_part,(data,host_address,offset,size))
     !>
     !> \ingroup grp_queues
     SUBROUTINE hmpprt_enqueue_upload_part(queue,data,host_address,offset,size)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T),    INTENT(IN) :: queue !< ...
       TYPE(HMPPRT_DATA_T),     INTENT(IN) :: data !< ...
       CLASS(*),                INTENT(IN) :: host_address(*) !< ...
       INTEGER(HMPPRT_SIZE_T),  INTENT(IN) :: offset   !< ...
       INTEGER(HMPPRT_SIZE_T),  INTENT(IN) :: size   !< ...
     END SUBROUTINE hmpprt_enqueue_upload_part
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
DUMMY_INTERFACE
     SUBROUTINE hmpprt_enqueue_upload_part(queue,data,host_address,offset,size)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T),    INTENT(IN) :: queue
       TYPE(HMPPRT_DATA_T),     INTENT(IN) :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  host_address
       INTEGER,                 INTENT(IN) :: host_address(*)
       INTEGER(HMPPRT_SIZE_T),  INTENT(IN) :: offset
       INTEGER(HMPPRT_SIZE_T),  INTENT(IN) :: size
     END SUBROUTINE hmpprt_enqueue_upload_part
  END_DUMMY_INTERFACE
#else
  EXTERNAL hmpprt_enqueue_upload_part
#endif


#ifdef DOC_MODE
  DUMMY_INTERFACE
     !>
     !> Enqueue the execution of REFCALL(hmpprt_data_download,(data,host_address))
     !>
     !> \ingroup grp_queues
    SUBROUTINE hmpprt_enqueue_download(queue,data,host_address)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN)  :: queue           !< ...
       TYPE(HMPPRT_DATA_T),  INTENT(IN)  :: data            !< ...
       CLASS(*),             INTENT(OUT) :: host_address(*) !< ...
     END SUBROUTINE hmpprt_enqueue_download
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE
    SUBROUTINE hmpprt_enqueue_download(queue,data,host_address)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN)  :: queue
       TYPE(HMPPRT_DATA_T),  INTENT(IN)  :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  host_address
       INTEGER,              INTENT(OUT) :: host_address(*)
     END SUBROUTINE hmpprt_enqueue_download
  END_DUMMY_INTERFACE
#else
  EXTERNAL hmpprt_enqueue_download
#endif

  !>
  !> Enqueue the execution of REFCALL(hmpprt_data_download_part,(data,host_address,offset,size))
  !> \ingroup grp_queues
  !>
#ifdef DOC_MODE
  DUMMY_INTERFACE
     SUBROUTINE hmpprt_enqueue_download_part(queue,data,host_address,offset,size)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN)       :: queue           !< ...
       TYPE(HMPPRT_DATA_T),  INTENT(IN)       :: data            !< ...
       CLASS(*),                INTENT(INOUT) :: host_address(*) !< ...
       INTEGER(HMPPRT_SIZE_T),  INTENT(IN)    :: offset          !< ...
       INTEGER(HMPPRT_SIZE_T),  INTENT()      :: size            !< ...
     END SUBROUTINE hmpprt_enqueue_download_part
  END_DUMMY_INTERFACE
#elif defined(DEC_NO_ARG_CHECK_SUPPORT)
  DUMMY_INTERFACE
     SUBROUTINE hmpprt_enqueue_download_part(queue,data,host_address,offset,size)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T),    INTENT(IN)    :: queue
       TYPE(HMPPRT_DATA_T),     INTENT(IN)    :: data
       !DEC$ ATTRIBUTES NO_ARG_CHECK ::  host_address
       INTEGER,                 INTENT(INOUT) :: host_address(*)
       INTEGER(HMPPRT_SIZE_T),  INTENT(IN)    :: offset
       INTEGER(HMPPRT_SIZE_T),  INTENT(IN)    :: size
     END SUBROUTINE hmpprt_enqueue_download_part
  END_DUMMY_INTERFACE
#else
  EXTERNAL hmpprt_enqueue_download_part
#endif

  DUMMY_INTERFACE
     !>
     !>  Enqueue the execution of REFCALL(hmpprt_queue_execute,(subqueue))
     !> \ingroup grp_queues
     !>
     SUBROUTINE hmpprt_enqueue_execute(queue,subqueue)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue    !< ...
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: subqueue !< ...
     END SUBROUTINE hmpprt_enqueue_execute
  END_DUMMY_INTERFACE


  DUMMY_INTERFACE
     !>
     !>  Enqueue the execution of REFCALL(hmpprt_queue_start,(subqueue))
     !> \ingroup grp_queues
     !>
     SUBROUTINE hmpprt_enqueue_start(queue,subqueue)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue !< ...
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: subqueue !< ...
     END SUBROUTINE hmpprt_enqueue_start
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !>  Enqueue the execution of REFCALL(hmpprt_queue_wait,(subqueue))
     !> \ingroup grp_queues
     !>
     SUBROUTINE hmpprt_enqueue_wait(queue,subqueue)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue    !< ...
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: subqueue !< ...
     END SUBROUTINE hmpprt_enqueue_wait
  END_DUMMY_INTERFACE

  DUMMY_INTERFACE
     !>
     !> Enqueue the execution of REFCALL(hmpprt_queue_clear,(subqueue))
     !> \ingroup grp_queues
     !>
     SUBROUTINE hmpprt_enqueue_clear(queue,subqueue)
       USE hmpprt_enum
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: queue    !< ...
       TYPE(HMPPRT_QUEUE_T), INTENT(IN) :: subqueue !< ...
     END SUBROUTINE hmpprt_enqueue_clear
  END_DUMMY_INTERFACE

END MODULE hmpprt
