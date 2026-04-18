#pragma once
#include "LinkedList.h"
#include "GenomicPosition.h"
#include "BaseLUT.h"
#include <string>
#include <cstdint>

namespace dna{

//Gere les LinkedList (buckets) et represente la table de dispersion   
class HashTable {
public:
    HashTable(const int& p_capacity);
    ~HashTable();
    
    bool encodeKmer(const std::string& p_kmer, uint64_t& p_encoded) const;
    uint8_t encodeBase(char p_base) const;
    void insert(const uint64_t& p_encoded, const GenomicPosition& p_pos);
    LinkedList* lookup(const uint64_t& p_encoded) const;  
    
private:
    LinkedList** m_buckets;
    int m_capacity;
};

}
