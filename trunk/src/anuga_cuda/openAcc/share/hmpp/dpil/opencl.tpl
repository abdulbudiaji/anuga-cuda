HEADER = (
#ifndef __GPUCODE__
#include \<stdio.h>
#include \<stdlib.h>
#include \<math.h>

#include \<hmpprt/Grouplet.h>
#include \<hmpprt/HostTypes.h>
#include \<hmpprt/Context.h>
#include \<hmpprt/OpenCLTypes.h>
#include \<hmpprt/OpenCLGrid.h>
#include \<hmpprt/OpenCLModule.h>
#include \<hmpprt/DeviceManager.h>
#include \<hmpperr/hmpperr.h>

#include \<CL/cl.h>

#ifdef _WIN32
#  define CDLT_API __declspec(dllexport)
#else /* ! _WIN32 */
#  define CDLT_API
#endif /* _WIN32 */

#else

#if defined(CL_KHR_FP64_SUPPORTED) && (defined(CL_VERSION_1_0) || defined(CL_VERSION_1_1))
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#endif

#ifdef GLOBAL_ATOMIC_EXTS_SUPPORTED
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics: enable
#endif

#ifdef LOCAL_ATOMIC_EXTS_SUPPORTED
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics: enable
#endif

#ifdef BYTE_ADDRESSABLE_STORE_EXTS_SUPPORTED
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store: enable
#endif

#endif /* __GPUCODE__ */

<forward_declarations>

#ifndef __GPUCODE__
<host_type_definitions>
#else
<dev_type_definitions>
<include_paths>
#endif
)

DEVICE_FUNCTION = (
#ifdef __GPUCODE__
<type> <name>(<params separator=", ">)
)

DEVICE_FUNCTION_PROTO = (
#ifdef __GPUCODE__
<type> <name>(<params separator=", ">);
)


GRID_FUNCTION_PROTO = (
#ifndef __GPUCODE__
static hmpprt::OpenCLGrid * <name> = 0;
#else
__kernel <kernel_attributes> void <name>(<params separator=", ">);
#endif // __GPUCODE__
)


GRID_FUNCTION = (
#ifdef __GPUCODE__
<const_params separator="\n">
__kernel <kernel_attributes> void <name>(<params separator=", ">)
)

AFTER_GRID_FUNCTION = #endif /* __GPUCODE__ */\n\n\n\n

FUNCTION = (
#ifndef __GPUCODE__
<linkage end=" ">CDLT_API <type> <name>(<params separator=", ">)
)

AFTER_FUNCTION = #endif /* __GPUCODE__ */\n\n\n\n

INIT_GRID_FUNCTION = <name> = new hmpprt::OpenCLGrid(hmpprt_module, "<name>");\n

FINI_GRID_FUNCTION = delete <name>;\n

HMPPRT_INIT = (
#ifndef __GPUCODE__
extern "C" const char * hmpprt_opencl_get_gpu_code();

static hmpprt::OpenCLModule * hmpprt_module = 0;
static int hmpprt_uses = 0;

extern "C" CDLT_API void * hmpprt_init()
{
  try
  {
    if (hmpprt_uses++ == 0)
    {
      hmpprt_module = new hmpprt::OpenCLModule(hmpprt_opencl_get_gpu_code());
      <first_time_content>
    }
    <every_time_content>
  }
  catch (hmpperr::Error & e)
  {
    return e.clone();
  }
  catch(...)
  {
    fprintf(stderr,"Unexpected error in hmpprt_init()\\n");
    abort();
  }
  return 0;
}
#endif /* __GPUCODE__ */
)

HMPPRT_FINI = (
#ifndef __GPUCODE__
extern "C" CDLT_API void * hmpprt_fini()
{
  try
  {
    if (--hmpprt_uses == 0)
    {
      <content>
      delete hmpprt_module;
      hmpprt_module = 0;
    }
  }
  catch (hmpperr::Error & e)
  {
    return e.clone();
  }
  catch(...)
  {
    fprintf(stderr,"Unexpected error in hmpprt_fini()\\n");
    abort();
  }
  return 0;
}
#endif /* __GPUCODE__ */
)

SHARED_MEM_DECLARE = (
local <type> <name><shape>;
)
