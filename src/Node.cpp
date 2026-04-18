#include "Node.h"

using namespace std;

namespace dna{
    Node::Node(const uint64_t& p_kmer, const GenomicPosition& p_pos):
        m_kmerCode(p_kmer), m_position(p_pos), m_next(nullptr) {};
}
