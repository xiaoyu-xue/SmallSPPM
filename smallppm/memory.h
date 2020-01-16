#pragma once
#define HAVE_ALIGNAS
#define L1_CACHE_LINE_SIZE 64
#include "def.h"
#include <list>
#include <algorithm>
#include "scalar.h"

NAMESPACE_BEGIN

#define ARENA_ALLOC(arena, Type) new ((arena).Alloc(sizeof(Type))) Type


void* AllocAligned(size_t size);

template <typename T>
T* AllocAligned(size_t count) {
    return (T*)AllocAligned(count * sizeof(T));
}

void FreeAligned(void*);


class
#ifdef HAVE_ALIGNAS
    alignas(L1_CACHE_LINE_SIZE)
#endif
    MemoryArena {
public:
    MemoryArena(size_t blockSize = 262144) : blockSize(blockSize) {}

    ~MemoryArena() {
        FreeAligned(currentBlock);
        for (auto& block : usedBlocks) FreeAligned(block.second);
        for (auto& block : availableBlocks) FreeAligned(block.second);
    }

    void* Alloc(size_t nBytes) {
        // Round up _nBytes_ to minimum machine alignment

#if !defined(HAVE_ALIGNOF)
        const int align = 16;
#else
        const int align = alignof(std::max_align_t);
#endif

        assert(IsPowerOf2(align));

        nBytes = (nBytes + align - 1) & ~(align - 1);
        if (currentBlockPos + nBytes > currentAllocSize) {
            // Add current block to _usedBlocks_ list
            if (currentBlock) {
                usedBlocks.push_back(
                    std::make_pair(currentAllocSize, currentBlock));
                currentBlock = nullptr;
                currentAllocSize = 0;
            }

            // Get new block of memory for _MemoryArena_

            // Try to get memory block from _availableBlocks_
            for (auto iter = availableBlocks.begin();
                iter != availableBlocks.end(); ++iter) {
                if (iter->first >= nBytes) {
                    currentAllocSize = iter->first;
                    currentBlock = iter->second;
                    availableBlocks.erase(iter);
                    break;
                }
            }
            if (!currentBlock) {
                currentAllocSize = std::max(nBytes, blockSize);
                currentBlock = AllocAligned<uint8_t>(currentAllocSize);
            }
            currentBlockPos = 0;
        }
        void* ret = currentBlock + currentBlockPos;
        currentBlockPos += nBytes;
        return ret;
    }
private:
    MemoryArena(const MemoryArena&) = delete;
    MemoryArena& operator=(const MemoryArena&) = delete;
    // MemoryArena Private Data
    const size_t blockSize;
    size_t currentBlockPos = 0, currentAllocSize = 0;
    uint8_t* currentBlock = nullptr;
    std::list<std::pair<size_t, uint8_t*>> usedBlocks, availableBlocks;
};

NAMESPACE_END