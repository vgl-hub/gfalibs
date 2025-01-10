#ifndef STRUCT
#define STRUCT

struct UserInput { // a container for user input

    // memory
    uint64_t maxMem = 0;
    // files
    std::string inSequence; // input file to evaluate
    std::string inSak; // input of instructions for the swiss army knife
    std::string inAgp; // input agp
    std::string inBedInclude; // input bed file of coordinates to include
    std::string inBedExclude; // input bed file of coordinates to exclude
    std::vector<std::string> inFiles; // input reads to evaluate
    std::string inAlign;
    // coordinates
    BedCoordinates bedIncludeList;
    // options
    char pipeType = 'n'; // default pipe type null
    std::string sortType = "none"; // type of sorting (default: none)
    
    int noSequence = 0, // output gfa without sequence
        discoverPaths_flag = 0,
        stats_flag = 0;
    
    std::string file(char type, uint16_t fileNum = 0);
    uint8_t mode = 0;
    int8_t hc_cutoff = -1;

    std::vector<std::string> kmerDB; // a database of kmers (or DBG)
    std::string prefix = ".", outFile = "";
    uint8_t kmerLen = 21;
    
    uint64_t gSize = 0; // expected genome size, with 0 statistics are not computed
    
    std::vector<std::string> outFiles; // output files
    uint32_t splitLength = 0;
    
};

struct Sequence { // a generic sequence container
    std::string header, comment;
    std::string sequence, sequenceQuality;
    unsigned int seqPos = 0;
};

struct Sequences { // a collection of sequences
    std::vector<Sequence> sequences;
    uint32_t batchN;
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

enum variantType {REF, SNV, INS, DEL, COM};
struct DBGpath {
    
    variantType type;
    uint16_t refLen = 1;
    uint64_t pos;
    std::string sequence;
    double score = 0;
    
};

#endif //STRUCT
