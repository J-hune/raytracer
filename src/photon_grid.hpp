#pragma once

#include <cuda_runtime.h>

#include <cstddef>

namespace rt {

struct Photon;

// Scratch bytes CUB needs to scan the per-cell counters.
std::size_t photonScanStorage(unsigned int bucketCount);

// Sorts the photons emitted by the tracer into contiguous per-cell runs. On
// return offsets[bucket] .. offsets[bucket + 1] delimits the photons of a cell
// inside sorted, and counts has been clobbered.
cudaError_t buildPhotonGrid(const Photon* photons, unsigned int count, float cellSize,
                            unsigned int bucketCount, unsigned int* counts,
                            unsigned int* offsets, void* scan, std::size_t scanBytes,
                            Photon* sorted);

}
