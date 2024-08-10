#ifndef STRUCT
#define STRUCT

#include "memory.h"

struct UserInput { // a container for user input

    // memory
    uint64_t maxMem = 0;
    // files
    std::string inSequence; // input file to evaluate
    std::string inSak; // input of instructions for the swiss army knife
    std::string inAgp; // input agp
    std::string inBedInclude; // input bed file of coordinates to include
    std::string inBedExclude; // input bed file of coordinates to exclude
    std::vector<std::string> inReads; // input reads to evaluate
    std::string inAlign;
    // coordinates
    BedCoordinates bedIncludeList;
    // options
    char pipeType = 'n'; // default pipe type null
    std::string sortType = "none"; // type of sorting (default: none)
    
    int noSequence = 0; // output gfa without sequence
    
    std::string file(char type, uint16_t fileNum = 0);
    uint8_t mode = 0,
            hc_cutoff = -1,
            discoverPaths_flag = 0,
            stats_flag = 0;
    
    std::vector<std::string> kmerDB; // a database of kmers (or DBG)
    std::string prefix = ".", outFile = "";
    
    uint32_t kmerLen = 21;
    
};

struct Sequence { // a generic sequence container
    
    std::string header, comment;
    std::string* sequence = NULL, *sequenceQuality = NULL;
    unsigned int seqPos = 0;
    
    ~Sequence();
    
};

struct Sequences { // a collection of sequences
    
    std::vector<Sequence*> sequences;
    unsigned int batchN;
    
    ~Sequences();
    
};

struct Tag {
    
    char type, label[3] = "";
    std::string content;
    
};

struct Gap {
    
    char orientation0;
    unsigned int segmentId;
    char orientation1;
    uint64_t dist;
    unsigned int edgeId;
    
};

struct Edge {
    
    char orientation0;
    unsigned int id;
    char orientation1;
    unsigned int weight = 0;
    
    bool operator==(const Edge& e) const;
    
};

enum ComponentType {SEGMENT, GAP, EDGE};
struct PathComponent {
    
    ComponentType componentType;
    uint32_t id;
    char orientation;
    uint64_t start;
    uint64_t end;
    
};

struct Bubble {
    unsigned int id0, id1, id2, id3;
};

enum variantType {REF, SNV, INS, DEL};
struct DBGpath {
    
    variantType type;
    uint64_t pos;
    std::string sequence;
    double score = 0;
    
    DBGpath() {}
    DBGpath(uint64_t pos) : pos(pos) {}
    DBGpath(uint64_t pos, double score) : pos(pos), score(score) {}
    DBGpath(variantType type, uint64_t pos, std::string sequence, double score) : type(type), pos(pos), sequence(sequence), score(score) {}
    
};

template<typename TYPE> // this is a generic buffer, TYPE is the type of the elements we wish to store in it. Usually each hashed kmer becomes part of a buffer specified by its hash value
struct Buf {
    uint64_t pos = 0, size; // pos keeps track of the position reached filling the buffer, initialized to contain up to size elements
    TYPE *seq = new TYPE[size](); // the actual container, parentheses ensure bits are set to 0
    
    Buf() : size(pow(2,8)){
        alloc += size*sizeof(TYPE);
    }
    Buf(uint64_t size) : size(size){
        alloc += size*sizeof(TYPE);
    }
    ~Buf(){
        if (seq != NULL) {
            delete[] seq;
            freed += size*sizeof(TYPE);
        }
    }
    
    uint64_t newPos(uint64_t add) {
        
        if (pos + add > size) {
            
            uint64_t newSize = (pos + add)*2;
            
            alloc += newSize*sizeof(TYPE);
            TYPE* seqNew = new TYPE[newSize](); // parentheses ensure bits are set to 0
            
            memcpy(seqNew, seq, size*sizeof(TYPE));
            
            delete[] seq;
            freed += size*sizeof(TYPE);
            size = newSize;
            seq = seqNew;
            
        }
        return pos += add;
    }
};

struct Buf2bit {
    uint64_t size;
    uint8_t *seq = new uint8_t[size](); // the actual container, parentheses ensure bits are set to 0
    
    Buf2bit(uint64_t size) : size(size){
        alloc += size*sizeof(uint8_t);
    }
    ~Buf2bit(){
        if (seq != NULL) {
            delete[] seq;
            freed += size*sizeof(uint8_t);
        }
    }
    inline uint8_t getBase(uint64_t index) const {
        uint8_t offset = index%4;
        return (seq[index/4] >> (6-offset*2)) & 3;
    }
};

#endif //STRUCT
