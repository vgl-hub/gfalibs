#include <iostream>
#include <vector>

struct GraphNode {
    
    uint8_t base;
    std::vector<GraphNode*> next;
    
    GraphNode(uint8_t base) : base(base) {}
    
    void addNext(GraphNode* nextNode) {
        next.push_back(nextNode);
    }
    
    bool exists(GraphNode* candidateNode) {
        for (GraphNode *nextNode : next)
            if (candidateNode->base == nextNode->base)
                return true;
        return false;
    }
    
};

struct StringGraph {
    
    GraphNode *root = new GraphNode(5);
    uint8_t *seq;
    uint8_t k;
    uint64_t pos;
    
    StringGraph(uint8_t *seq, uint8_t k, uint64_t start = 0) : seq(seq), k(k), pos(start) {
        
        GraphNode *prev = root;
        
        for (uint8_t i = 0; i < k; ++i) {
            GraphNode *graphNode = new GraphNode(seq[pos]);
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
            prev = graphNode;
        }
        pos += k;
        
    }

    uint64_t currentPos() {
        
        return pos;
        
    }
    
    void advancePos(uint64_t n) {
        
        pos += n;
        
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
    
    void appendNext(uint8_t n = 1) {
        
        for (uint8_t i = 0; i < n; ++i) {
            GraphNode* newLeaf = new GraphNode(seq[pos++]);
            addNext(root, newLeaf);
        }
    }
    
    void addNext(GraphNode* nextNode, GraphNode* newLeaf) {
        
        if (nextNode->next.size() != 0) {
            for (GraphNode *nextGraphNode : nextNode->next)
                addNext(nextGraphNode, newLeaf);
        }else{
            if (nextNode != newLeaf && !nextNode->exists(newLeaf))
                nextNode->addNext(newLeaf);
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
                if (!nextNode->exists(newLeaf))
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
                for (uint8_t i = 0; i < k + 2; ++i) // so that we can check k+2 kmers downstream if needed
                    newPath.push_back(seq[pos+i]);
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
