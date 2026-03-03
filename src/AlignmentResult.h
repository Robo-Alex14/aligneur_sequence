#pragma once
#include "GenomicPosition.h"

// Result of aligning a read to one genomic window.
class AlignmentResult {
public:
    GenomicPosition position;
    int             mismatches;   // positions that differ (excluding non-ACGT)
    int             matches;      // positions that agree  (excluding non-ACGT)

    AlignmentResult();
    AlignmentResult(const GenomicPosition& pos, int mismatches, int matches);

    bool operator<(const AlignmentResult& o) const;
    bool operator>(const AlignmentResult& o) const;
};
