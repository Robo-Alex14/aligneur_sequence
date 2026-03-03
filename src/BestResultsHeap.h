#pragma once
#include "AlignmentResult.h"
#include <cstdint>

// Fixed-capacity max-heap ordered by mismatch count.
//
// We keep the N BEST (lowest-mismatch) results seen so far, where N is the
// capacity passed to the constructor (default 15).
// A max-heap: the root is always the WORST result
// AMONG those currently stored, so we can check in O(1)
// whether a new candidate is worth
// keeping, and evict the root in O(log N) when it is.
//
// 
// The underlying storage is a raw heap-allocated array of AlignmentResult.
//
class BestResultsHeap {
public:
    //a vous de jouer
private:
    //a vous de jouer
};
