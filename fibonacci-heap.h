//
//  fibonacci-heap.h
//  kreeq-dev
//
//  modified from https://arxiv.org/abs/2303.10034
//

#ifndef fibonacci_heap_h
#define fibonacci_heap_h

//Struct used for each Fibonacci heap node
template<typename V>
struct FibonacciNode {
    int degree;
    FibonacciNode<V>* parent;
    FibonacciNode<V>* child;
    FibonacciNode<V>* left;
    FibonacciNode<V>* right;
    bool mark;
    int key;
    V objPtr;
};
//Fibonacci heap class
template<typename V>
class FibonacciHeap {
    FibonacciNode<V>* minNode;
    uint16_t numNodes, maxNumNodes = 1000;
    std::vector<FibonacciNode<V>*> degTable;
    phmap::parallel_flat_hash_map<uint64_t, FibonacciNode<V>*> nodePtrs; // this originally was a vector
    public:
    FibonacciHeap() {
        //Constructor function
        this->numNodes = 0;
        this->minNode = NULL;
        this->degTable = {};
    }
    ~FibonacciHeap() {
        //Destructor function
        this->numNodes = 0;
        this->minNode = NULL;
        for (auto node : nodePtrs)
            delete node.second;
        this->degTable.clear();
        this->nodePtrs.clear();
    }
    int size() {
        //Number of nodes in the heap
        return this->numNodes;
    }
    bool empty() {
        //Is the heap empty?
        if (this->numNodes > 0) return false;
        else return true;
    }
    void insert(V u, int key) {
        if (numNodes >= maxNumNodes) { // make space for the new node
            decreaseKey(degTable.back()->objPtr, 0);
            V toBeDeleted = extractMin();
            nodePtrs.erase(toBeDeleted->first);
            delete toBeDeleted;
        }
        //Insert the vertex u with the specified key (value for L(u)) into the Fibonacci heap. O(1) operation
        this->nodePtrs[u->first] = new FibonacciNode<V>;
        this->nodePtrs[u->first]->objPtr = u;
        FibonacciNode<V>* node = this->nodePtrs[u->first];
        node->key = key;
        node->degree = 0;
        node->parent = NULL;
        node->child = NULL;
        node->left = node;
        node->right = node;
        node->mark = false;
        FibonacciNode<V>* minN = this->minNode;
        if (minN != NULL) {
            FibonacciNode<V>* minLeft = minN->left;
            minN->left = node;
            node->right = minN;
            node->left = minLeft;
            minLeft->right = node;
        }
        if (minN == NULL || minN->key > node->key) {
            this->minNode = node;
        }
        this->numNodes++;
    }
    V extractMin() {
        //Extract the node with the minimum key from the heap. O(log n) operation, where n is the number of nodes in the heap
        FibonacciNode<V>* minN = this->minNode;
        if (minN != NULL) {
//            std::cout<<"here we are"<<std::endl;
            int deg = minN->degree;
            FibonacciNode<V>* currChild = minN->child;
            FibonacciNode<V>* remChild;
            for (int i = 0; i < deg; i++) {
//                std::cout<<"here we are0.1"<<std::endl;
                remChild = currChild;
//                std::cout<<"here we are0.2"<<std::endl;
                currChild = currChild->right;
                _existingToRoot(remChild);
//                std::cout<<"here we are0.3"<<std::endl;
            }
//            std::cout<<"here we are1"<<std::endl;
            _removeNodeFromRoot(minN);
//            std::cout<<"here we are2"<<std::endl;
            this->numNodes--;
            if (this->numNodes == 0) {
                this->minNode = NULL;
            }else{
//                std::cout<<"here we are2.1"<<std::endl;
                this->minNode = minN->right;
                FibonacciNode<V>* minNLeft = minN->left;
                this->minNode->left = minNLeft;
//                std::cout<<"here we are2.2"<<std::endl;
                minNLeft->right = this->minNode;
                _consolidate();
//                std::cout<<"here we are2.3"<<std::endl;
            }
            return minN->objPtr;
            
        }else{
//            std::cout<<"here we are3"<<std::endl;
            return NULL;
        }
        
    }
    void decreaseKey(V u, int newKey) {
        //Decrease the key of the node in the Fibonacci heap that has index u. O(1) operation
        FibonacciNode<V>* node = this->nodePtrs[u->first];
        if (newKey > node->key) return;
        node->key = newKey;
        if (node->parent != NULL) {
            if (node->key < node->parent->key) {
                FibonacciNode<V>* parentNode = node->parent;
                _cut(node);
                _cascadingCut(parentNode);
            }
        }
        if (node->key < this->minNode->key) {
            this->minNode = node;
        }
    }
    private:
    //The following are private functions used by the public methods above
    void _existingToRoot(FibonacciNode<V>* newNode) {
        FibonacciNode<V>* minN = this->minNode;
        newNode->parent = NULL;
        newNode->mark = false;
        if (minN != NULL) {
            FibonacciNode<V>* minLeft = minN->left;
            minN->left = newNode;
            newNode->right = minN;
            newNode->left = minLeft;
            minLeft->right = newNode;
            if (minN->key > newNode->key) {
                this->minNode = newNode;
            }
        }
        else {
            this->minNode = newNode;
            newNode->right = newNode;
            newNode->left = newNode;
        }
    }
    void _removeNodeFromRoot(FibonacciNode<V>* node) {
        if (node->right != node) {
            node->right->left = node->left;
            node->left->right = node->right;
        }
        if (node->parent != NULL) {
            if (node->parent->degree == 1) {
                node->parent->child = NULL;
            }
            else {
                node->parent->child = node->right;
            }
            node->parent->degree--;
        }
    }
    void _cut(FibonacciNode<V>* node) {
        _removeNodeFromRoot(node);
        _existingToRoot(node);
    }
    void _addChild(FibonacciNode<V>* parentNode, FibonacciNode<V>* newChildNode) {
        if (parentNode->degree == 0) {
            parentNode->child = newChildNode;
            newChildNode->right = newChildNode;
            newChildNode->left = newChildNode;
            newChildNode->parent = parentNode;
        }
        else {
            FibonacciNode<V>* child1 = parentNode->child;
            FibonacciNode<V>* child1Left = child1->left;
            child1->left = newChildNode;
            newChildNode->right = child1;
            newChildNode->left = child1Left;
            child1Left->right = newChildNode;
        }
        newChildNode->parent = parentNode;
        parentNode->degree++;
    }
    void _cascadingCut(FibonacciNode<V>* node) {
        FibonacciNode<V>* parentNode = node->parent;
        if (parentNode != NULL) {
            if (node->mark == false) {
                node->mark = true;
            }
            else {
                _cut(node);
                _cascadingCut(parentNode);
            }
        }
    }
    void _link(FibonacciNode<V>* highNode, FibonacciNode<V>* lowNode) {
        _removeNodeFromRoot(highNode);
        _addChild(lowNode, highNode);
        highNode->mark = false;
    }
    void _consolidate() {
        int deg, rootCnt = 0;
//        std::cout<<"hey1"<<std::endl;
        if (this->numNodes > 1) {
            this->degTable.clear();
            FibonacciNode<V>* currNode = this->minNode;
            FibonacciNode<V>* currDeg, * currConsolNode;
            FibonacciNode<V>* temp = this->minNode, * itNode = this->minNode;
//            std::cout<<"hey2"<<std::endl;
            do {
                rootCnt++;
                itNode = itNode->right;
            } while (itNode != temp);
//            std::cout<<"hey3"<<std::endl;
            for (int cnt = 0; cnt < rootCnt; cnt++) {
                currConsolNode = currNode;
                currNode = currNode->right;
                deg = currConsolNode->degree;
//                std::cout<<"hey4"<<std::endl;
                while (true) {
//                    std::cout<<"hey4.1"<<std::endl;
                    while (deg >= int(this->degTable.size())) {
                        this->degTable.push_back(NULL);
                    }
                    if (this->degTable[deg] == NULL) {
                        this->degTable[deg] = currConsolNode;
                        break;
                    }else{
//                        std::cout<<"hey4.2"<<std::endl;
                        currDeg = this->degTable[deg];
                        if (currConsolNode->key > currDeg->key) {
                            std::swap(currConsolNode, currDeg);
                        }
                        if (currDeg == currConsolNode) break;
                        _link(currDeg, currConsolNode);
                        this->degTable[deg] = NULL;
                        deg++;
//                        std::cout<<"hey4.3"<<std::endl;
                    }
//                    std::cout<<"hey5"<<std::endl;
                }
//                std::cout<<"hey6"<<std::endl;
            }
            this->minNode = NULL;
//            std::cout<<"hey7"<<std::endl;
            for (size_t i = 0; i < this->degTable.size(); i++) {
                if (this->degTable[i] != NULL) {
                    _existingToRoot(this->degTable[i]);
                }
            }
//            std::cout<<"hey8"<<std::endl;
        }
    }
};

#endif /* fibonacci_heap_h */
