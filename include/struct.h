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
    
    uint32_t kLen = 21;
    uint8_t sLen = 0;
    
    uint64_t gSize = 0; // expected genome size, with 0 statistics are not computed
    
    std::vector<std::string> outFiles; // output files
    uint32_t splitLength = 0;
    
};

struct Sequence { // a generic sequence container
    
    std::string header, comment;
    std::string* sequence = NULL, *sequenceQuality = NULL;
    unsigned int seqPos = 0;
    
    ~Sequence();
    
    void deleteSequence();
    
};

struct Sequences { // a collection of sequences
    
    std::vector<Sequence*> sequences;
    uint32_t batchN;
    
    ~Sequences();
    
};

struct Sequence2 {
	std::string header, comment;
	std::unique_ptr<std::string> sequence;        // replace raw
	std::unique_ptr<std::string> sequenceQuality; // replace raw
	unsigned int seqPos = 0;
	// deleteSequence() no longer needed
	Sequence2() {}
	Sequence2(std::string h, std::string c, std::unique_ptr<std::string> s)
	  : header(std::move(h)), comment(std::move(c)), sequence(std::move(s)) {}
	Sequence2(std::string h, std::string c, std::unique_ptr<std::string> s, std::unique_ptr<std::string> q)
	  : header(std::move(h)), comment(std::move(c)), sequence(std::move(s)), sequenceQuality(std::move(q)) {}
	
	void set(std::string h, std::string c,
			 std::unique_ptr<std::string> s,
			 std::unique_ptr<std::string> q = nullptr,
			 unsigned int p = 0) {
		header  = std::move(h);
		comment = std::move(c);
		sequence = std::move(s);
		sequenceQuality = std::move(q);
		seqPos = p;
	}
};

struct Sequences2 {
	std::vector<std::unique_ptr<Sequence2>> sequences; // replace raw
	uint32_t batchN = 0, fileN = 0;
	
	void recycle_keep_capacity() {
		for (auto& up : sequences) {
			if (!up) continue;                 // slot is empty (nullptr)
			Sequence2& s = *up;
			s.header.clear();                  // keep capacity
			s.comment.clear();
			if (s.sequence)        s.sequence->clear();         // keep capacity
			if (s.sequenceQuality) s.sequenceQuality->clear();  // keep capacity
			s.seqPos = 0;
		}
		sequences.resize(0);
	}
	
	void add(uint32_t seqPos, uint32_t batchCounter, std::string seqHeader, std::string seqComment, std::unique_ptr<std::string> inSequence, std::unique_ptr<std::string> inSequenceQuality = std::unique_ptr<std::string>()) {
		if (sequences.size() <= batchCounter)
			sequences.emplace_back(std::make_unique<Sequence2>());

		sequences[batchCounter]->set(
			std::move(seqHeader),
			std::move(seqComment),
			std::move(inSequence),
			std::move(inSequenceQuality),
			seqPos
		);
	}
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

template<typename TYPE = uint8_t> // this is a generic buffer, TYPE is the type of the elements we wish to store in it. Usually each hashed kmer becomes part of a buffer specified by its hash value
struct Buf {
    uint64_t pos = 0, size = 0; // pos keeps track of the position reached filling the buffer, initialized to contain up to size elements
    TYPE *seq = NULL;
    
    Buf() : size(pow(2,8)){
        alloc += size*sizeof(TYPE);
        seq = new TYPE[size](); // the actual container, parentheses ensure bits are set to 0
    }
    Buf(uint64_t size) : size(size){
        alloc += size*sizeof(TYPE);
        seq = new TYPE[size](); // the actual container, parentheses ensure bits are set to 0
    }
    Buf(const Buf& buf) {
        pos = buf.pos;
        size = buf.size;
        seq = buf.seq;
    }
    Buf(Buf&& buf) {
        pos = std::move(buf.pos);
        size = std::move(buf.size);
        seq = std::move(buf.seq);
        buf.seq = NULL;
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

#endif //STRUCT
