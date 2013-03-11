HEADER = (
#include \<stdio.h>

#ifndef __CUDACC__
#include \<stdlib.h>
#include \<math.h>

#include \<hmpprt/Grouplet.h>
#include \<hmpprt/HostTypes.h>
#include \<hmpprt/Context.h>
#include \<hmpprt/CUDAGrid.h>
#include \<hmpprt/CUDAModule.h>
#include \<hmpprt/DeviceManager.h>
#include \<hmpperr/hmpperr.h>

#ifdef _WIN32
#  define CDLT_API __declspec(dllexport)
#else /* ! _WIN32 */
#  define CDLT_API
#endif /* _WIN32 */

<dwarf_xml>

#else // ! __CUDACC__

#include \<hmpprt/HostTypes.h>
#include \<hmpprt/CUDAIntrinsics.h>

extern __shared__ int64_t hmpp_sharedmem[];
#endif // __CUDACC__

<forward_declarations>

#ifndef __CUDACC__
<host_type_definitions>
#else
<dev_type_definitions>
<include_paths>
#endif

)


DEVICE_FUNCTION = (
#ifdef __CUDACC__
<linkage end=" ">__device__ <type> <name>(<params separator=", ">)
)

GRID_FUNCTION_PROTO = (
#ifndef __CUDACC__
static hmpprt::CUDAGrid * <name> = 0;
#else
<const_params separator="\n">
extern "C" __global__ void <name>(<params separator=", ">);
#endif // __CUDACC__\n\n\n\n
)

GRID_FUNCTION = (
#ifdef __CUDACC__
<const_params separator="\n">
extern "C" __global__ void <name>(<params separator=", ">)
)

AFTER_GRID_FUNCTION = #endif // __CUDACC__\n\n\n\n

FUNCTION = (
#ifndef __CUDACC__
<linkage end=" "><type> <name>(<params separator=", ">)
)

AFTER_FUNCTION = #endif // __CUDACC__\n\n\n\n

INIT_GRID_FUNCTION = <name> = new hmpprt::CUDAGrid(hmpprt_module, "<name>");\n

FINI_GRID_FUNCTION = delete <name>;\n

HMPPRT_INIT = (
#ifndef __CUDACC__
extern "C" const char * hmpprt_cuda_get_gpu_code();

static hmpprt::CUDAModule * hmpprt_module = 0;
static int hmpprt_uses = 0;

extern "C" CDLT_API void * hmpprt_init()
{
  try
  {
    if (hmpprt_uses++ == 0)
    {
      hmpprt_module = new hmpprt::CUDAModule(hmpprt_cuda_get_gpu_code());
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
#endif // __CUDACC__
)

HMPPRT_FINI = (
#ifndef __CUDACC__
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
#endif // __CUDACC__
)

SHARED_MEM_DECLARE = (
<type> * <name> = (<type> *)(((char *)hmpp_sharedmem + <offset>));
)

CONST_MEM_INNER_DECLARE = (
typedef <type> <name>_t <shape>;
<name>_t & <name> = (<name>_t &)(<funcname>_const_<name>);
)


CONSTADDR_MEM_INNER_DECLARE = (
<type> * <name> = <funcname>_const_<name>;
)


CONST_MEM_OUTER_DECLARE = __constant__ <type> <stars> <funcname>_const_<name> <shape>;

