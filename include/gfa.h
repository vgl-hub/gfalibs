//
//  gfastats-gfa.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFA_H
#define GFA_H

class InSequences { //collection of InSegment and inGap objects and their summary statistics
    
protected:

    std::vector<Log> logs;
    
    //gfa variables
    std::vector<InSegment*> inSegments;
    std::vector<InGap> inGaps;
    std::vector<InEdge> inEdges;
    std::vector<InPath> inPaths;
    std::vector<std::vector<Gap>> adjListFW;
    std::vector<std::vector<Gap>> adjListBW;
    std::vector<std::vector<Edge>> adjEdgeList;
    phmap::flat_hash_map<std::string, unsigned int> headersToIds;
    phmap::flat_hash_map<unsigned int, std::string> idsToHeaders;
    phmap::flat_hash_map<int, bool> visited, deleted;
    bool backward = false;
    
    std::vector<uint64_t> scaffLens;
    std::vector<uint64_t> contigLens;
    std::vector<uint64_t> gapLens;
    
    std::vector<uint64_t> scaffNstars   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffLstars   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<uint64_t> scaffNGstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffLGstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<uint64_t> contigNstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigLstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<uint64_t> contigNGstars {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigLGstars {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<uint64_t> gapNstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> gapLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    double scaffAuN = 0, scaffAuNG = 0, contigAuN = 0, contigAuNG = 0, gapAuN = 0;
    
    InSegment inSegment;
    InGap gap;
    InEdge edge;
    InPath path;
    
    uint64_t
    totScaffLen = 0,
    totContigLen = 0,
    totSegmentLen = 0,
    totGapLen = 0,
    totA = 0,
    totC = 0,
    totG = 0,
    totT = 0,
    totLowerCount = 0;
    
    unsigned int
    scaffN = 0,
    contigN = 0,
    pathN = 0;

    //connectivity
    unsigned int deadEnds = 0;
    unsigned int disconnectedComponents = 0;
    uint64_t lengthDisconnectedComponents = 0;
    
    std::vector<Bubble> bubbles;
    
    friend class SAK;
    friend class Threadpool;
    
public:
    
    ~InSequences();
    
    UIdGenerator uId; // unique numeric identifier for each feature
    
    std::vector<Log> getLogs();
    
    InGap pushbackGap(Log* threadLog, InPath* path, std::string* seqHeader, unsigned int* iId, unsigned int* dist, char sign, unsigned int uId1, unsigned int uId2);
    
    InSegment* pushbackSegment(unsigned int currId, Log* threadLog, InPath* path, std::string* seqHeader, std::string* seqComment, std::string* sequence, unsigned int* iId, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint64_t sStart, uint64_t sEnd, std::string* sequenceQuality = NULL);
    
    bool traverseInSequence(Sequence* sequence, int hc_cutoff);
    
    bool traverseInSegment(Sequence* sequence, std::vector<Tag> inSequenceTags);
    
    void appendSequence(Sequence* sequence, int hc_cutoff);
    
    void appendSegment(Sequence* sequence, std::vector<Tag> inSequenceTags);
    
    InSegment *getInSegment(unsigned int sId);
    
    std::vector<InSegment*>* getInSegments();
    
    std::vector<InGap>* getInGaps();
    
    InPath getInPath(unsigned int pId);
    
    std::vector<InPath> getInPaths();
    
    uint64_t getTotScaffLen();
    
    uint64_t getTotSegmentLen();
    
    void changeTotGapLen(unsigned int gapLen);
    
    unsigned int getGapNScaffold();
    
    unsigned int getGapN();

    unsigned int getEdgeN();
    
    unsigned int getPathN();
    
    void recordScaffLen(uint64_t seqLen);
    
    void recordGapLen(unsigned int gapLen);
    
    void evalNstars(char type, uint64_t gSize = 0);
    
    void evalAuN(char type, uint64_t gSize = 0);
    
    void computeAuN(std::vector<uint64_t>& lens, double& auN, double* auNG = 0, uint64_t gSize = 0);
    
    unsigned int getScaffN();
    
    unsigned int getTotContigN();
    
    std::vector <uint64_t> getScaffNstars();
    
    std::vector <uint64_t> getScaffNGstars();
    
    std::vector <unsigned int> getScaffLstars();
    
    std::vector <unsigned int> getScaffLGstars();
    
    std::vector <uint64_t> getContigNstars();
    
    std::vector <uint64_t> getContigNGstars();
    
    std::vector <unsigned int> getContigLstars();
    
    std::vector <unsigned int> getContigLGstars();
    
    std::vector <uint64_t> getGapNstars();
    
    std::vector <unsigned int> getGapLstars();
    
    uint64_t getScaffN50();
    
    uint64_t getScaffNG50();
    
    unsigned int getScaffL50();
    
    unsigned int getScaffLG50();
    
    double getScaffauN();
    
    double getScaffauNG();
    
    double getContigauN();
    
    double getContigauNG();
    
    double getGapauN();
    
    unsigned int getSegmentN();
    
    unsigned int getContigN50();
    
    unsigned int getContigNG50();
    
    unsigned int getContigL50();
    
    unsigned int getContigLG50();
    
    unsigned int getGapN50();
    
    unsigned int getGapL50();
    
    uint64_t getLargestScaffold();

    uint64_t getSmallestScaffold();

    uint64_t getLargestContig();

    uint64_t getSmallestContig();

    unsigned int getLargestGap();

    unsigned int getSmallestGap();
    
    double computeAvgScaffLen();
    
    uint64_t getTotContigLen();

    double computeAvgContigLen();
    
    double computeAvgSegmentLen();
    
    uint64_t getTotGapLen();
    
    double computeAverageGapLen();
    
    uint64_t getTotA();
    
    uint64_t getTotC();
    
    uint64_t getTotG();
    
    uint64_t getTotT();
    
    uint64_t getTotLowerCount();
    
    double computeGCcontent();
    
    //gfa methods
    bool addGap(InGap inGap);
    
    bool addPath(InPath path);
    
    std::vector<InGap> getGaps();

    std::vector<InEdge>* getEdges();
    
    bool appendEdge(InEdge edge);
    
    //sorting methods

    void sortSegmentsByOriginal();
    
    void sortPathsByOriginal();
    
    void sortPathsByNameAscending();
    
    void sortPathsByNameDescending();
    
    void sortPathsByList(std::vector<std::string> headerList);
    
    void sortPathsBySize(bool largest);
    
    //gfa methods
    void insertHash(const std::string &segHeader, unsigned int i);
    
    unsigned int getuId();
    
    phmap::flat_hash_map<std::string, unsigned int>* getHash1();

    phmap::flat_hash_map<unsigned int, std::string>* getHash2();
    
    void buildGraph(std::vector<InGap> const& gaps);

    void buildEdgeGraph(std::vector<InEdge> const& edges);

    void dfsEdges(unsigned int v, unsigned int* componentLength);
    
    void dfsScaffolds(unsigned int v, unsigned int* scaffSize, unsigned int* A, unsigned int* C, unsigned int* G, unsigned int* T, unsigned int* lowerCount); // Depth First Search to explore graph connectivity
    
    std::vector<std::vector<Gap>> getAdjListFW();
    
    std::vector<std::vector<Gap>> getAdjListBW();
    
    bool getVisited(unsigned int uId);
    
    bool getDeleted(unsigned int uId);
    
    bool updateStats();
    
    bool removeTerminalGaps();

    unsigned int getDeadEnds();

    unsigned int getDisconnectedComponents();
    
    unsigned int getLengthDisconnectedComponents();
    
    // instruction methods
    
    std::vector<InGap> getGap(std::string* contig1, std::string* contig2 = NULL);
    
    std::vector<unsigned int> removeGaps(std::string* contig1, std::string* contig2 = NULL);
    
    bool deleteSegment(std::string* contig1);
    
    void removePath(unsigned int pUId, bool all = false, bool silent = false);
    
    void removeGap(unsigned int gUId, bool silent = false);
    
    void resizeGap(std::string gHeader, unsigned int size);
    
    void removePathsFromSegment(unsigned int uId);
    
    void removePathComponents(unsigned int uId);
    
    void removeSegmentInPath(unsigned int suId, InGap gap);
    
    void joinPaths(std::string pHeader, unsigned int pUId1, unsigned int pUId2, std::string gHeader, unsigned int gUId, char pId1Or, char pId2Or, unsigned int dist, unsigned int start1, unsigned int end1, unsigned int start2, unsigned int end2);
    
    InPath joinPathsByComponent(std::string seqHeader, unsigned int uId1, unsigned int uId2, unsigned int uId3);
    
    void splitPath(unsigned int guId, std::string pHeader1, std::string pHeader2);
    
    void clearPaths();
    
    void clearGaps();
    
    void renamePath(unsigned int pUId, std::string pHeader, unsigned int* newpUId);
    
    void revComPath(unsigned int pUId);
    
    void trimPathByUId(unsigned int pUId, unsigned int start, unsigned int end);

    void trimPathByRef(std::vector<PathComponent>& pathComponents, unsigned int start, unsigned int end);

    void trimPath(std::vector<PathComponent>* pathComponents, unsigned int start, unsigned int end);
    
    void trimComponent(PathComponent& component, int start, int end);
    
    int getComponentSize(PathComponent& component, bool original);
    
    unsigned int pathLen(unsigned int pUId);
    
    void walkPath(InPath* path);
    
    void discoverPaths();
    
    unsigned int detectOverlap(std::string* sequence1, std::string* sequence2, unsigned int minOvlLen);
    
    void discoverTerminalOverlaps(int terminalOvlLen);
    
    void dfsPath(unsigned int v, InPath& newPath); // Depth First Search to build a new path given a vertex
    
    void findBubbles();
    
    std::vector<Bubble>* getBubbles();
    
    std::vector<unsigned int> getCircular();
    
    void maskPath(std::string pHeader, unsigned int start, unsigned int end, unsigned int dist);
    
};

#endif /* GFA_H */
