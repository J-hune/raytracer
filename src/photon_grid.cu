#include "photon_grid.hpp"

#include "gpu_shared.hpp"
#include "photon_hash.cuh"

#include <cub/device/device_scan.cuh>

namespace rt {
namespace {

constexpr unsigned int BlockSize = 256U;

unsigned int blockCount(unsigned int count) {
    return (count + BlockSize - 1U) / BlockSize;
}

__global__ void countPhotons(const Photon* photons, unsigned int count, float cellSize,
                             unsigned int bucketCount, unsigned int* counts) {
    const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count)
        return;
    const unsigned int bucket =
        photonBucket(photonCell(photons[index].position, cellSize), bucketCount);
    atomicAdd(&counts[bucket], 1U);
}

__global__ void scatterPhotons(const Photon* photons, unsigned int count, float cellSize,
                               unsigned int bucketCount, const unsigned int* offsets,
                               unsigned int* cursors, Photon* sorted) {
    const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count)
        return;
    const Photon photon = photons[index];
    const unsigned int bucket =
        photonBucket(photonCell(photon.position, cellSize), bucketCount);
    sorted[offsets[bucket] + atomicAdd(&cursors[bucket], 1U)] = photon;
}

}

std::size_t photonScanStorage(unsigned int bucketCount) {
    std::size_t bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, bytes, static_cast<unsigned int*>(nullptr),
                                  static_cast<unsigned int*>(nullptr), bucketCount + 1U);
    return bytes;
}

cudaError_t buildPhotonGrid(const Photon* photons, unsigned int count, float cellSize,
                            unsigned int bucketCount, unsigned int* counts,
                            unsigned int* offsets, void* scan, std::size_t scanBytes,
                            Photon* sorted) {
    const std::size_t counterBytes = (bucketCount + 1U) * sizeof(unsigned int);
    cudaError_t status = cudaMemset(counts, 0, counterBytes);
    if (status != cudaSuccess)
        return status;
    if (count) {
        countPhotons<<<blockCount(count), BlockSize>>>(photons, count, cellSize,
                                                       bucketCount, counts);
        if ((status = cudaGetLastError()) != cudaSuccess)
            return status;
    }

    // Turns the per-cell counts into the start offset of each cell's run. The
    // extra trailing element makes offsets[bucket + 1] a valid end for the last
    // bucket.
    status = cub::DeviceScan::ExclusiveSum(scan, scanBytes, counts, offsets,
                                           bucketCount + 1U);
    if (status != cudaSuccess)
        return status;

    // The counts have been consumed by the scan, so they double as the write
    // cursors of the scatter below.
    if ((status = cudaMemset(counts, 0, counterBytes)) != cudaSuccess)
        return status;
    if (count) {
        scatterPhotons<<<blockCount(count), BlockSize>>>(photons, count, cellSize,
                                                         bucketCount, offsets, counts,
                                                         sorted);
        if ((status = cudaGetLastError()) != cudaSuccess)
            return status;
    }
    return cudaDeviceSynchronize();
}

}
