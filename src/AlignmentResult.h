#pragma once
#include "GenomicPosition.h"

namespace dna {

class AlignmentResult {
public:
    GenomicPosition m_position;
    int             m_mismatches;   
    int             m_matches;      

    AlignmentResult();
    AlignmentResult(const GenomicPosition& pos, int mismatches, int matches);

    bool operator<(const AlignmentResult& o) const;
    bool operator>(const AlignmentResult& o) const;
};

}
