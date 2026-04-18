#include "LinkedList.h"

namespace dna{
    
    LinkedList::LinkedList() : m_head(nullptr){};
    
    //Elimination de tout les noeuds contenus dans la linked list
    LinkedList::~LinkedList() {
        Node* current = m_head;
        while(current != nullptr){
            Node* temp = current;
            current = current->m_next;
            delete temp;
        }
    }
    
    void LinkedList::insert(const uint64_t& p_kmerCode, const GenomicPosition& p_pos){
        Node* insertion = new Node(p_kmerCode, p_pos);
        insertion->m_next = m_head;
        m_head = insertion;
    };
    
    Node* LinkedList::first() const{
        return m_head;
    };
    
}
