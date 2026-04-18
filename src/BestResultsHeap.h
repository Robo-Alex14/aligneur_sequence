#pragma once
#include "AlignmentResult.h"
#include <cstdint>
#include <vector>

namespace dna {
    
class BestResultsHeap {
public:
    BestResultsHeap(const int& p_capacity = 15);

    void insert(const AlignmentResult& result);
    const AlignmentResult getBest() const;

    bool isEmpty() const;
    bool isFull() const;
    int getCapacity() const;
    size_t getSize() const;
    
private:
    
    std::vector<AlignmentResult> m_data;
    int m_capacity;
    
    //Methodes internes
    void siftUp(int p_index);
    void siftDown(int p_index);
};

}
