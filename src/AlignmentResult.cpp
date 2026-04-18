#include "AlignmentResult.h"

namespace dna {

AlignmentResult::AlignmentResult()
    : m_position(), m_mismatches(0), m_matches(0) {};

AlignmentResult::AlignmentResult(const GenomicPosition& p_pos, int p_mismatches, 
        int p_matches)
    : m_position(p_pos), m_mismatches(p_mismatches), m_matches(p_matches) {};

bool AlignmentResult::operator<(const AlignmentResult& o) const {
    return m_mismatches < o.m_mismatches;
};

bool AlignmentResult::operator>(const AlignmentResult& o) const {
    return m_mismatches > o.m_mismatches;
};
}
