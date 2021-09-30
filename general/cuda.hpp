// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_CUDA_HPP
#define MFEM_CUDA_HPP

#include "../config/config.hpp"
#include "error.hpp"
#include "device.hpp"
#include <iostream>

// CUDA block size used by MFEM.
#define MFEM_CUDA_BLOCKS 256

#ifdef MFEM_USE_CUDA
#define MFEM_DEVICE __device__
#define MFEM_LAMBDA __host__
#define MFEM_HOST_DEVICE __host__ __device__
#define MFEM_DEVICE_SYNC MFEM_GPU_CHECK(cudaDeviceSynchronize())
#define MFEM_STREAM_SYNC MFEM_GPU_CHECK(cudaStreamSynchronize(0))
#define  MFEM_TIME_CALL(cnt_,msg_,...) CuTimeCall<cnt_>(msg_,[&]() {__VA_ARGS__})

// Define a CUDA error check macro, MFEM_GPU_CHECK(x), where x returns/is of
// type 'cudaError_t'. This macro evaluates 'x' and raises an error if the
// result is not cudaSuccess.
#define MFEM_GPU_CHECK(x) \
   do \
   { \
      cudaError_t err = (x); \
      if (err != cudaSuccess) \
      { \
         mfem_cuda_error(err, #x, _MFEM_FUNC_NAME, __FILE__, __LINE__); \
      } \
   } \
   while (0)
#endif // MFEM_USE_CUDA

// Define the MFEM inner threading macros
#if defined(MFEM_USE_CUDA) && defined(__CUDA_ARCH__)
#define MFEM_SHARED __shared__
#define MFEM_SYNC_THREAD __syncthreads()
#define MFEM_BLOCK_ID(k) blockIdx.k
#define MFEM_THREAD_ID(k) threadIdx.k
#define MFEM_THREAD_SIZE(k) blockDim.k
#define MFEM_FOREACH_THREAD(i,k,N) for(int i=threadIdx.k; i<N; i+=blockDim.k)
#endif

namespace mfem
{

#ifdef MFEM_USE_CUDA
// Function used by the macro MFEM_GPU_CHECK.
void mfem_cuda_error(cudaError_t err, const char *expr, const char *func,
                     const char *file, int line);
#endif

/// Allocates device memory and returns destination ptr.
void* CuMemAlloc(void **d_ptr, size_t bytes);

/// Allocates managed device memory
void* CuMallocManaged(void **d_ptr, size_t bytes);

/// Allocates page-locked (pinned) host memory
void* CuMemAllocHostPinned(void **ptr, size_t bytes);

/// Frees device memory and returns destination ptr.
void* CuMemFree(void *d_ptr);

/// Frees page-locked (pinned) host memory and returns destination ptr.
void* CuMemFreeHostPinned(void *ptr);

/// Copies memory from Host to Device and returns destination ptr.
void* CuMemcpyHtoD(void *d_dst, const void *h_src, size_t bytes);

/// Copies memory from Host to Device and returns destination ptr.
void* CuMemcpyHtoDAsync(void *d_dst, const void *h_src, size_t bytes);

/// Copies memory from Device to Device
void* CuMemcpyDtoD(void *d_dst, const void *d_src, size_t bytes);

/// Copies memory from Device to Device
void* CuMemcpyDtoDAsync(void *d_dst, const void *d_src, size_t bytes);

/// Copies memory from Device to Host
void* CuMemcpyDtoH(void *h_dst, const void *d_src, size_t bytes);

/// Copies memory from Device to Host
void* CuMemcpyDtoHAsync(void *h_dst, const void *d_src, size_t bytes);

/// Check the error code returned by cudaGetLastError(), aborting on error.
void CuCheckLastError();

/// Get the number of CUDA devices
int CuGetDeviceCount();

template <int count, typename LMBDA>
inline void CuTimeCall(const char msg[], LMBDA &&lmbda)
{
   cudaEvent_t start,stop;
   float       timeT = 0.0f,timeLo = 9999999999.9f,timeHi = 0.0f;
   int         runCount = 0;

   MFEM_GPU_CHECK(cudaEventCreate(&start));
   MFEM_GPU_CHECK(cudaEventCreate(&stop));
   MFEM_DEVICE_SYNC;
   for (; runCount < count; ++runCount)
   {
      float time = 0.0f;

      MFEM_GPU_CHECK(cudaEventRecord(start,0));
      lmbda();
      MFEM_GPU_CHECK(cudaEventRecord(stop,0));
      MFEM_GPU_CHECK(cudaEventSynchronize(stop));
      MFEM_GPU_CHECK(cudaEventElapsedTime(&time,start,stop));
      timeLo = timeLo <= time ? timeLo : time;
      timeHi = timeHi >= time ? timeHi : time;
      timeT += time;
   }
   MFEM_GPU_CHECK(cudaEventDestroy(start));
   MFEM_GPU_CHECK(cudaEventDestroy(stop));
   std::cout<<msg<<"\n";
   std::cout<<"Total number of iterations "<<runCount<<"\n";
   std::cout<<"Total elapsed time "<<timeT<<" (ms)\n";
   std::cout<<"Average time/run "<<timeT/(double)std::max(runCount,
                                                          1)<<" (ms) Lo "<<timeLo<<" (ms), Hi "<<timeHi<<" (ms)\n";
   return;
}

template <std::size_t N> class DeviceStreamPool<Backend::CUDA,N>;

#ifdef MFEM_USE_CUDA
template <std::size_t N>
class DeviceStreamPool<Backend::CUDA,N>
{
   using stream_type = cudaStream_t;
   std::array<stream_type,N> _pool;

public:
   DeviceStreamPool()
   {
      for (auto&& stream : _pool) { MFEM_GPU_CHECK(cudaStreamCreate(&stream)); }
   }

   ~DeviceStreamPool()
   {
      for (auto&& stream : _pool) { MFEM_GPU_CHECK(cudaStreamDestroy(stream)); }
   }

   stream_type get(int idx = MFEM_STREAM_NONE) const
   {
      static int currentStream = 0;

      switch (idx)
      {
         case MFEM_STREAM_NONE:
            return nullptr;
         case MFEM_STREAM_NEXT:
            currentStream = ++currentStream % N;
            return _pool[currentStream];
         default:
            MFEM_ASSERT(idx >= 0 &&
                        idx < N,"idx "<<idx<<" does not satisfy 0 <=  "idx<<" < "<<N);
            return _pool[idx];
      }
      return nullptr;
   }
};

const DeviceStreamPool<Backend::CUDA,5>& cudaPool();
#endif

} // namespace mfem

#endif // MFEM_CUDA_HPP
