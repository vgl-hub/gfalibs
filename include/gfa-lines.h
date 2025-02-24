#ifndef GFA_LINES_H
#define GFA_LINES_H

class InSegment { // DNA sequence with no gaps
protected:
    std::string seqHeader;
    std::string seqComment;
    std::string* inSequence = NULL;
    std::string* inSequenceQuality = NULL;
    uint64_t A = 0, C = 0, G = 0, T = 0, N = 0, lowerCount = 0;
    unsigned int uId = 0, iId = 0, seqPos = 0;
    std::vector<Tag> tags;
    std::vector<std::deque<DBGpath>> variants;
    
    friend class SAK;
    friend class InSequences;
    friend class Report;
    
public:
    
    InSegment() {};
    
    InSegment(InSegment &inSegment) {
        seqHeader = inSegment.seqHeader;
        seqComment = inSegment.seqComment;
        if (inSegment.inSequence != NULL)
            inSequence = new std::string(*inSegment.inSequence);
        if (inSegment.inSequenceQuality != NULL)
            inSequenceQuality = new std::string(*inSegment.inSequenceQuality);
        A = inSegment.A, C = inSegment.C, G = inSegment.G, T = inSegment.T, N = inSegment.N, lowerCount = inSegment.lowerCount;
        uId = inSegment.uId, iId = inSegment.iId, seqPos = inSegment.seqPos;
        tags = inSegment.tags;
    };
    
    ~InSegment();
    
    void set(Log* threadLog, unsigned int uId, unsigned int iId, std::string seqHeader, std::string* seqComment, std::string* sequence, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, unsigned int seqPos, std::string* sequenceQuality = NULL, std::vector<Tag>* inSequenceTags = NULL, uint64_t* N = NULL);
    
    void setSeqHeader(std::string h);
    
    void setSeqComment(std::string c);
    
    void setInSequence(std::string* s);
    
    void setInSequenceQuality(std::string* q);
    
    void setSeqTags(std::vector<Tag>* t);

    void setuId(unsigned int i); // absolute id
    
    void setiId(unsigned int i); // temporary id, internal to scaffold
    
    void setSeqPos(unsigned int i);

    std::string getSeqHeader();
    
    std::string getSeqComment();
    
    std::vector<Tag> getTags();
    
    std::string getInSequence(unsigned int start = 0, unsigned int end = 0);
    
    std::string* getInSequencePtr();
    
    std::string getInSequenceQuality(unsigned int start = 0, unsigned int end = 0);
    
    unsigned int getSeqPos();
    
    uint64_t getSegmentLen(uint64_t start = 0, uint64_t end = 0);
    
    unsigned int getuId(); // absolute id
    
    unsigned int getsId(); // temporary id, internal to scaffold
    
    void setACGT(uint64_t* a, uint64_t* c, uint64_t* g, uint64_t* t, uint64_t* n = NULL);
    
    void setLowerCount(uint64_t* C);
    
    uint64_t getA();
    
    uint64_t getC();
    
    uint64_t getG();
    
    uint64_t getT();

    uint64_t getN();
    
    unsigned int getLowerCount(uint64_t start = 0, uint64_t end = 0);
    
    double computeGCcontent();
    
    bool trimSegment(uint64_t start, uint64_t end);
    
    bool updateSegmentCounts(uint64_t start, uint64_t end);
    
    bool rvcpSegment();
    
    bool invertSegment();
    
    bool isCircularSegment(std::vector<unsigned int>* circularSegments);
    
    char* first();
    
    void addVariants(std::vector<std::deque<DBGpath>> variants);
    
    std::vector<std::deque<DBGpath>>& getVariants();
    
};

class InGap {
private:
//    uint64_t lineN; // useful if we wish to sort as is the original input
    std::string gHeader;
    char sId1Or, sId2Or;
    unsigned int uId, iId, sId1, sId2;
    uint64_t dist;
    std::vector<Tag> tags;
    
    friend class SAK;
    friend class InSequences;
    
public:
    void newGap(unsigned int uId, unsigned int sId1, unsigned int sId2, const char& sId1or, const char& sId2or, uint64_t dist, std::string gHeader = "", std::vector<Tag> tags = {});

    void setuId(unsigned int i); // absolute id
    
    void setiId(unsigned int i); // temporary id, internal to scaffold

    void setsId1(unsigned int i);
    
    void setsId2(unsigned int i);
    
    void setDist(uint64_t i);
    
    std::string getgHeader();
    
    unsigned int getuId();
    
    unsigned int getsId1();
    
    char getsId1Or();
    
    unsigned int getsId2();
    
    char getsId2Or();
    
    uint64_t getDist(uint64_t start = 0, uint64_t end = 0);
    
    std::vector<Tag> getTags();
    
    
    
};
class InEdge {
    private:
//    uint64_t lineN; // useful if we wish to sort as is the original input
    unsigned int eUId = 0, eId, sId1, sId2;
    char sId1Or, sId2Or;
    std::string cigar, eHeader;
    std::vector<Tag> tags;
    
    friend class SAK;
    friend class InSequences;
    
public:
    
    InEdge() {};
 
    InEdge(unsigned int eUId, unsigned int eId, unsigned int sId1, unsigned int sId2, const char& sId1Or, const char& sId2Or, std::string cigar = "", std::string eHeader = "", std::vector<Tag> tags = {}) : eUId(eUId), eId(eId), sId1(sId1), sId2(sId2), sId1Or(sId1Or), sId2Or(sId2Or), cigar(cigar), eHeader(eHeader), tags(tags) {};
    
    InEdge(unsigned int sId1, unsigned int sId2, const char& sId1Or, const char& sId2Or, std::string cigar = "", std::string eHeader = "", std::vector<Tag> tags = {}) : sId1(sId1), sId2(sId2), sId1Or(sId1Or), sId2Or(sId2Or), cigar(cigar), eHeader(eHeader), tags(tags) {};
    
    void newEdge(unsigned int eUId, unsigned int sId1, unsigned int sId2, const char& sId1Or, const char& sId2Or, std::string cigar = "", std::string eHeader = "", std::vector<Tag> tags = {});
    
    bool operator==(const InEdge& e) const;

    void seteUId(unsigned int i); // absolute id
    
    void seteId(unsigned int i); // temporary id, internal to scaffold

    void setsId1(unsigned int i);
    
    void setsId2(unsigned int i);
    
    void setSeqTags(std::vector<Tag>* t);
    
    void appendTag(Tag t);
    
    std::string geteHeader();
    
    std::string getCigar();
    
    unsigned int geteUId();

    unsigned int geteId();
    
    unsigned int getsId1();
    
    char getsId1Or();
    
    unsigned int getsId2();
    
    char getsId2Or();
    
    std::vector<Tag> getTags();
    
    bool isSelf();
    
};

class InPath {
    
private:
//    uint64_t lineN; // useful if we wish to sort as is the original input
    std::string pHeader, pComment;
    std::vector<PathComponent> pathComponents;
    uint32_t pUId, contigN = 0, seqPos;
    
    uint64_t length = 0, segmentLength = 0, lowerCount = 0, A = 0, C = 0, G = 0, T = 0;
    
    friend class SAK;
    friend class InSequences;

public:
    
    void newPath(unsigned int pUid, std::string h, std::string c = "", unsigned int seqpos = 0);

    void setpUId(unsigned int pUid);
    
    void setHeader(std::string pheader);
    
    void setComment(std::string c);
    
    void add(ComponentType componentType, unsigned int UId, char sign = '+', uint64_t start = 0, uint64_t end = 0);
    
    void append(std::vector<PathComponent> components);
    
    void insert(std::vector<PathComponent>::iterator prevComponent, PathComponent newComponent);
    
    void clearPath();
    
    void setComponents(std::vector<PathComponent> newComponents);
    
    std::vector<PathComponent> getComponents();
    
    std::vector<PathComponent>* getComponentsByRef();

    unsigned int getpUId();
    
    std::string getHeader();
    
    std::string getComment();
    
    unsigned int getSeqPos();
    
    unsigned int getContigN();
    
    uint64_t getLen();
    
    uint64_t getA();
    
    uint64_t getC();
    
    uint64_t getG();
    
    uint64_t getT();
    
    uint64_t getSegmentLen();
    
    uint64_t getLowerCount();
    
    void revCom();
    
    void increaseContigN();
    
    void increaseGapN();
    
    void increaseLen(uint64_t n);
    
    void increaseSegmentLen(uint64_t n);
    
    void increaseLowerCount(uint64_t n);
    
    void increaseA(uint64_t n);
    
    void increaseC(uint64_t n);
    
    void increaseG(uint64_t n);
    
    void increaseT(uint64_t n);
    
    void reinitializeCounts();
    
    bool isCircularPath(std::vector<uint32_t>* circularPaths);
    
};

#endif /* GFA_LINES_H */
