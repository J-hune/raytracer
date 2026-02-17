#pragma once

#include <cuda_runtime.h>

// Shared by the OptiX photon tracer and the CUDA grid builder: both have to
// agree on the cell decomposition, so the hash lives in one place.

namespace rt {

__device__ __forceinline__ int3 photonCell(float3 position, float cellSize) {
    return make_int3(
        static_cast<int>(floorf(position.x / cellSize)),
        static_cast<int>(floorf(position.y / cellSize)),
        static_cast<int>(floorf(position.z / cellSize)));
}

__device__ __forceinline__ unsigned int photonBucket(int3 cell,
                                                     unsigned int bucketCount) {
    const unsigned int hash =
        static_cast<unsigned int>(cell.x) * 73856093U ^
        static_cast<unsigned int>(cell.y) * 19349663U ^
        static_cast<unsigned int>(cell.z) * 83492791U;
    return hash & (bucketCount - 1U);
}

}
