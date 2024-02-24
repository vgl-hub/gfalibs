#include <iostream>
#include <vector>

struct GraphNode {
    
    uint8_t base;
    std::vector<GraphNode*> next;
    
    GraphNode(uint8_t base) : base(base) {}
    
    void addNext(GraphNode* nextNode) {
        next.push_back(nextNode);
    }
    
};

struct StringGraph {
    
    GraphNode *root = new GraphNode(5);
    std::vector<GraphNode*> leaves;
    uint8_t *seq;
    uint8_t k;
    uint64_t pos = k-1;
    GraphNode *prev;
    
    
    StringGraph(uint8_t *seq, uint8_t k) : seq(seq), k(k) {
        
        GraphNode *prev = root;
        
        for (uint8_t pos = 0; pos < k-1; ++pos) {
            GraphNode *graphNode = new GraphNode(seq[pos]);
            prev->addNext(graphNode);
            this->prev = prev;
            prev = graphNode;
        }
        leaves.push_back(prev);

    }

    uint8_t peek() {
        
        return seq[pos];
        
    }
    
    void pop_back() {
        
        if (root->next.size() > 0) {
            std::vector<GraphNode*> tmp = root->next[0]->next;
            for (GraphNode* nextNode : root->next)
                delete nextNode;
            root->next = tmp;
        }
    }
    
    void appendNode() {
        
        GraphNode* newLeaf = new GraphNode(seq[pos++]);
        for (GraphNode* leaf : leaves)
            leaf->addNext(newLeaf);
        leaves.clear();
        leaves.push_back(newLeaf);
        this->prev = newLeaf;
            
    }
    
    void addAlt(std::vector<uint8_t> alts) {
        
        leaves.clear();
        
        for (uint8_t alt : alts) {
            
            GraphNode *newLeaf = new GraphNode(alt);
            this->prev->addNext(newLeaf);
            leaves.push_back(newLeaf);
            
        }
        
        ++pos;
        pop_back();
            
    }
    
    void addIns(uint8_t altBase) {
        

    }
    
    std::vector<std::vector<uint8_t>> walkStringGraph(GraphNode *graphNode, std::vector<uint8_t> currentPath){
        
        std::vector<std::vector<uint8_t>> paths;
        if (graphNode->base != 4 && graphNode->base != 5 )
            currentPath.push_back(graphNode->base);
        
        for (GraphNode *nextGraphNode : graphNode->next) {
            if (nextGraphNode->next.size() == 0) {
                currentPath.push_back(nextGraphNode->base);
                currentPath.push_back(seq[pos]);
                currentPath.push_back(seq[pos+1]);
                paths.push_back(currentPath);
                break;
            }
            std::vector<std::vector<uint8_t>> discoveredPaths = walkStringGraph(nextGraphNode, currentPath);
            paths.insert(paths.end(), discoveredPaths.begin(), discoveredPaths.end());
        }
        return paths;
    }
    
    void deleteStringGraph(GraphNode *graphNode) {
        
        if (graphNode->next.size() > 0)
            deleteStringGraph(graphNode->next[0]);
        for (uint8_t n = 1; n < graphNode->next.size(); ++n)
            delete graphNode->next[n];
        
        delete graphNode;

    }
    
    ~StringGraph() {
        deleteStringGraph(root);
    }
    
};
