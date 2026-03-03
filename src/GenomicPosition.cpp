#include "GenomicPosition.h"

GenomicPosition::GenomicPosition()
    : chrom(""), offset(0) {}

GenomicPosition::GenomicPosition(const std::string& chrom, uint32_t offset)
    : chrom(chrom), offset(offset) {}
