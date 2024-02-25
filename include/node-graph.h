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
    uint8_t *seq;
    uint8_t k;
    uint64_t pos = 0;
    
    StringGraph(uint8_t *seq, uint8_t k) : seq(seq), k(k) {
        
        GraphNode *prev = root;
        
        for (uint8_t p = 0; p < k-1; ++p) {
            GraphNode *graphNode = new GraphNode(seq[p]);
            prev->addNext(graphNode);
            prev = graphNode;
            ++pos;
        }

    }
    
    void backtrack(uint8_t *seq, uint8_t k, uint64_t b){
        
        deleteStringGraph(root);
        
        root = new GraphNode(5);
        pos -= k+b;
        GraphNode *prev = root;
        
        for (uint8_t p = 0; p < k; ++p) {
            GraphNode *graphNode = new GraphNode(seq[pos+p]);
            prev->addNext(graphNode);
        }
        pos += k;
        
    }

    uint8_t currentPos() {
        
        return pos;
        
    }
    
    uint8_t peek() {
        
        return seq[pos];
        
    }
    
    void pop_front() {
        
        if (root->next.size() > 0) {
            std::vector<GraphNode*> tmp = root->next[0]->next;
            for (GraphNode* nextNode : root->next)
                delete nextNode;
            root->next = tmp;
        }
    }
    
    void appendNext() {
        
        GraphNode* newLeaf = new GraphNode(seq[pos++]);
        addNext(root, newLeaf);
        
    }
    
    void addNext(GraphNode* nextNode, GraphNode* newLeaf) {
        
//        std::cout<<+nextNode->base<<std::endl;
        
        if (nextNode->next.size() != 0) {
//            std::cout<<"size "<<nextNode->next.size()<<std::endl;
            for (GraphNode *nextGraphNode : nextNode->next)
                addNext(nextGraphNode, newLeaf);
        }else{
            if (nextNode != newLeaf)
                nextNode->addNext(newLeaf);
//            std::cout<<"end"<<std::endl;
        }
            
    }
    
    void appendAlts(std::vector<uint8_t> alts) {
        
        addAlt(root, alts);
        pop_front();
        ++pos;
            
    }
    
    void addAlt(GraphNode* nextNode, std::vector<uint8_t> alts) {
        
        if (nextNode->next.size() != 0) {
            addAlt(nextNode->next[0], alts);
        }else{
            for (uint8_t alt : alts) {
                GraphNode *newLeaf = new GraphNode(alt);
                nextNode->addNext(newLeaf);
            }
        }
    }
    
    std::vector<std::vector<uint8_t>> walkStringGraph(GraphNode *graphNode, std::vector<uint8_t> currentPath){
        
//        std::cout<<+graphNode->base<<std::endl;
        
        std::vector<std::vector<uint8_t>> paths;
        if (graphNode->base != 4 && graphNode->base != 5)
            currentPath.push_back(graphNode->base);
        
        for (GraphNode *nextGraphNode : graphNode->next) {
            if (nextGraphNode->next.size() == 0) {
                std::vector<uint8_t> newPath = currentPath;
                if (graphNode->base != 4)
                    newPath.push_back(nextGraphNode->base);
                newPath.push_back(seq[pos]);
                newPath.push_back(seq[pos+1]);
                newPath.push_back(seq[pos+2]);
                newPath.push_back(seq[pos+3]);
                paths.push_back(newPath);
            }else{
//                std::cout<<std::to_string(graphNode->next.size())<<std::endl;
                std::vector<std::vector<uint8_t>> discoveredPaths = walkStringGraph(nextGraphNode, currentPath);
                paths.insert(paths.end(), discoveredPaths.begin(), discoveredPaths.end());
            }
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
