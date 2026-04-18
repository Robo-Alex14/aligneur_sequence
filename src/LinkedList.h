#pragma once
#include "Node.h"
#include "GenomicPosition.h"
#include <cstdint>


namespace dna{

//Chaine de collision (contient/gerer la memoire des nodes) 
class LinkedList {
public:
    LinkedList();
    ~LinkedList();
    
    void insert(const uint64_t& p_kmer, const GenomicPosition& p_pos);
    Node* first() const;
    
private:    
    Node* m_head;
};

}
