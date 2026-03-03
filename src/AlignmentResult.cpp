#include "AlignmentResult.h"

AlignmentResult::AlignmentResult()
    : position(), mismatches(0), matches(0) {}

AlignmentResult::AlignmentResult(const GenomicPosition& pos, int mm, int m)
    : position(pos), mismatches(mm), matches(m) {}

bool AlignmentResult::operator<(const AlignmentResult& o) const {
    return mismatches < o.mismatches;
}
bool AlignmentResult::operator>(const AlignmentResult& o) const {
    return mismatches > o.mismatches;
}
