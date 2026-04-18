#include "BestResultsHeap.h"
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <vector>

using namespace std;

namespace dna {

BestResultsHeap::BestResultsHeap(const int& p_capacity)
    :m_capacity(p_capacity) {
    m_data.reserve(p_capacity);
}

bool BestResultsHeap::isEmpty() const {
    return m_data.empty();
}

bool BestResultsHeap::isFull() const {
    return static_cast<int>(m_data.size()) == m_capacity;
}

int BestResultsHeap::getCapacity() const {
    return m_capacity;
}

size_t BestResultsHeap::getSize() const{
    return m_data.size();
};

const AlignmentResult BestResultsHeap::getBest() const {
    if (isEmpty()) {
        throw std::runtime_error("Results are empty!");
    }

    size_t best_index = 0;
    for (size_t i = 1; i < m_data.size(); ++i) {
        if (m_data[i] < m_data[best_index]) {
            best_index = i;
        }
    }
    return m_data[best_index];
}

void BestResultsHeap::insert(const AlignmentResult& result) {
    if (!isFull()) {
        m_data.push_back(result);
        siftUp(m_data.size() - 1); //Get the last element added!
    }else if (result < m_data[0]) {
        m_data[0] = result;
        siftDown(0);
    }
}

void BestResultsHeap::siftUp(int p_index) {
    while (p_index > 0 && p_index < static_cast<int>(m_data.size())) {
        int parent = (p_index - 1) / 2;

        if (m_data[parent] < m_data[p_index]) {
            swap(m_data[parent], m_data[p_index]);
            p_index = parent;
        }else{
            break;
        }
    }
}

void BestResultsHeap::siftDown(int p_index) {
    int left = 2 * p_index + 1;
    int size = static_cast<int>(m_data.size());

    //Si index d'element a gauche est valide
    while (left < size) {
        int right = left + 1;
        int largest = p_index;

        if (m_data[largest] < m_data[left]) {
            largest = left;
        } //on verifie que right existe et qui ne depasse pas la limite de m_data
        else if (right < size && m_data[largest] < m_data[right]) {
            largest = right;
        }else if (largest == p_index) {
            break;
        }

        swap(m_data[p_index], m_data[largest]);
        p_index = largest;

        //calculer prochain element a gauche
        left = 2 * p_index + 1; 
    }
};

}
