#pragma once
#include <string>
#include <cstdint>

// Stores a chromosome/sequence name and a 0-based position within it.
class GenomicPosition {
public:
    std::string chrom;
    uint32_t    offset;

    GenomicPosition();
    GenomicPosition(const std::string& chrom, uint32_t offset);
};
