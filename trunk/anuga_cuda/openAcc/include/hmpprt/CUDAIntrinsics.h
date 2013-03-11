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
/// \internal

#ifndef HMPPRT_CUDA_INTRINSICS_H
#define HMPPRT_CUDA_INTRINSICS_H

#ifdef __CUDACC__

namespace hmpprt
{

__device__ __inline__ int gr_atidf()
{
  return (blockDim.x * blockDim.y * blockIdx.x  +  threadIdx.y * blockDim.x  +  threadIdx.x);
}

__device__ __inline__ int gr_atidx()
{
  return (blockDim.x * blockIdx.x  +  threadIdx.x);
}

__device__ __inline__ int gr_atidy()
{
  return (blockDim.y * blockIdx.y  +  threadIdx.y);
}

__device__ __inline__ int gr_btidf()
{
  return (threadIdx.y * blockDim.x  +  threadIdx.x);
}

__device__ __inline__ int gr_btidx()
{
  return threadIdx.x;
}

__device__ __inline__ int gr_btidy()
{
  return threadIdx.y;
}

__device__ __inline__ int gr_btnumf()
{
  return (blockDim.x * blockDim.y);
}

__device__ __inline__ int gr_btnumx()
{
  return blockDim.x;
}

__device__ __inline__ int gr_btnumy()
{
  return blockDim.y;
}

__device__ __inline__ int gr_gbidf()
{
  return (gridDim.x * blockIdx.y + blockIdx.x);
}

__device__ __inline__ int gr_gbidx()
{
  return blockIdx.x;
}

__device__ __inline__ int gr_gbidy()
{
  return blockIdx.y;
}

__device__ __inline__ int gr_gbnumf()
{
  return (gridDim.x * gridDim.y);
}

__device__ __inline__ int gr_gbnumx()
{
  return gridDim.x;
}

__device__ __inline__ int gr_gbnumy()
{
  return gridDim.y;
}

__device__ __inline__ int gr_atidminf()
{
  return blockDim.x * blockDim.y * blockIdx.x;
}

__device__ __inline__ int gr_atidminx()
{
  return blockDim.x * blockIdx.x;
}

__device__ __inline__ int gr_atidminy()
{
  return  blockDim.y * blockIdx.y;
}

__device__ __inline__ void gr_barrier()
{
  __syncthreads();
}

template <typename T>
__device__ __inline__ void gr_copyin(T * dst, T * src, size_t num_element)
{
  for (size_t i = 0 ; i < num_element ; i += gr_btnumf())
  {
    size_t pos = i + gr_btidf();
    if (pos < num_element)
    {
      dst[pos] = src[pos];
    }
  }
}

} // namespace hmpprt

#endif // __CUDACC__

#endif // HMPPRT_CUDA_INTRINSICS_H
