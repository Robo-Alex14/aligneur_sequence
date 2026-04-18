#include "HashTable.h"
#include "BaseLUT.h"

using namespace std;

namespace dna{
    
    //Peupler hashTable
    HashTable::HashTable(const int& p_capacity): m_capacity(p_capacity){
        m_buckets = new LinkedList*[m_capacity];
        
        for(int i = 0; i < m_capacity; i++){
            m_buckets[i] = new LinkedList();
        }
    };
    
    //Gestion de memoire - Eliminer les buckets
    HashTable::~HashTable(){
        for(int i = 0; i < m_capacity; i++){
            delete m_buckets[i];
        }
        delete[] m_buckets;
    };
    
    uint8_t HashTable::encodeBase(char p_base) const{
        return BASE_LUT[static_cast<unsigned char>(p_base)];
    };
    
    bool HashTable::encodeKmer(const std::string& p_kmer, uint64_t& p_encoded) const {
        p_encoded = 0;

        for (std::size_t i = 0; i < p_kmer.size(); ++i) {
            uint8_t modified_base = encodeBase(p_kmer[i]);

            if (modified_base == INVALID) {
                return false;
            }
            p_encoded = (p_encoded << 2) | modified_base;
        }
        return true;
    }
    
    void HashTable::insert(const uint64_t& p_encoded_kmer, const GenomicPosition& p_pos) {
        uint32_t index = static_cast<uint32_t>(p_encoded_kmer % m_capacity);
        m_buckets[index]->insert(p_encoded_kmer, p_pos); 
    }
    
    LinkedList* HashTable::lookup(const uint64_t& p_encoded_kmer) const{     
        uint32_t index = static_cast<uint32_t>(p_encoded_kmer % m_capacity);
        return m_buckets[index];
    };
}
