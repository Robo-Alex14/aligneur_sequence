#pragma once
#include "LinkedList.h"
#include "GenomicPosition.h"
#include "BaseLUT.h"
#include <string>
#include <cstdint>

// Hash table: kmer (2-bit encoded) -> linked list of (kmerCode, GenomicPosition).
// Bucket index  = kmerCode % capacity.
// Collision chains store the kmerCode in every Node so lookup() callers can
// skip nodes that belong to a different kmer sharing the same bucket.
//
// buckets_ is a raw heap-allocated array of LinkedList* (not a std::vector)
// so manage the allocation explicitly.
class HashTable {
public:
    // a vous de jouer
private:
    // a vous de jouer
};
