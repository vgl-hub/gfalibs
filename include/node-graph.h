#include <iostream>
#include <string>
#include <vector>

struct GraphNode {
    
    char base;
    std::vector<GraphNode*> next;
    
    GraphNode(char base) : base(base) {}
    
    void addNext(GraphNode* nextNode) {
        next.push_back(nextNode);
    }
    
};

void deleteStringGraph(GraphNode *graphNode) {
    
    if (graphNode->next.size() > 0)
        deleteStringGraph(graphNode->next[0]);
    delete graphNode;

}

std::vector<std::string> walkStringGraph(GraphNode *graphNode, std::string currentPath){
    
    std::vector<std::string> paths;
    if (graphNode->base != 'J' && graphNode->base != '0' )
        currentPath+=graphNode->base;
    
    for (GraphNode *nextGraphNode : graphNode->next) {
        if (nextGraphNode->next.size() == 0) {
            paths.push_back(currentPath+=nextGraphNode->base);
            break;
        }
        std::vector<std::string> discoveredPaths = walkStringGraph(nextGraphNode, currentPath);
        paths.insert(paths.end(), discoveredPaths.begin(), discoveredPaths.end());
    }
    return paths;
}

std::vector<std::string> stringToGraph(std::string seq, uint64_t err, uint8_t k) {
    
    GraphNode *root = new GraphNode('0');
    GraphNode *prev = root;
    
    for (uint64_t i = 0; i < seq.size()-k+1; ++i) {
        //std::cout<<seq.substr(i, k)<<std::endl;
        if (i+k-1 == err) {
            
            for (uint8_t pos = i+1; pos < i+k-1; ++pos) {
                GraphNode *graphNode = new GraphNode(seq[pos]);
                prev->addNext(graphNode);
                prev = graphNode;
            }

            GraphNode *pri = new GraphNode(seq[i+k-1]);
            prev->addNext(pri);
            GraphNode *alt = new GraphNode('A');
            prev->addNext(alt);
            GraphNode *alt2 = new GraphNode('J');
            prev->addNext(alt2);
            GraphNode *next = new GraphNode(seq[i+k]);
            pri->addNext(next);
            alt->addNext(next);
            alt2->addNext(next);
            prev = next;
            
        }
        
    }
    std::vector<std::string> paths = walkStringGraph(root, std::string());
    deleteStringGraph(root);
    return paths;
}
