#pragma once
#include "GenomicPosition.h"
#include <cstdint>

namespace dna{
    
//Represente un Kmer    
class Node {
    public:
        Node(const uint64_t& p_kmer, const GenomicPosition& p_pos);
     
        uint64_t m_kmerCode;
        GenomicPosition m_position;
        Node* m_next;
};
}

