#include <cstdint>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>

#include "log.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"

#include "gfa-lines.h"

InSegment::~InSegment()
{
    delete inSequence;
    delete inSequenceQuality;
}

void InSegment::set(Log* threadLog, unsigned int uId, unsigned int iId, std::string seqHeader, std::string* seqComment, std::string* sequence, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, unsigned int seqPos, std::string* sequenceQuality, std::vector<Tag>* inSequenceTags, uint64_t* N) {
    
    threadLog->add("Processing segment: " + seqHeader + " (uId: " + std::to_string(uId) + ", iId: " + std::to_string(iId) + ")");
    
    uint64_t seqSize = 0;
    
    this->setiId(iId); // set temporary sId internal to scaffold
    
    this->setuId(uId); // set absolute id
    
    this->setSeqPos(seqPos); // set original order
    
    this->setSeqHeader(seqHeader);
    
    if (*seqComment != "") {
        
        this->setSeqComment(*seqComment);
        
    }
    
    if (inSequenceTags != NULL) {
        
        this->setSeqTags(inSequenceTags);
        
    }
    
    if (*sequence != "*") {
        
        this->setInSequence(sequence);
        
        threadLog->add("Segment sequence set");
        
        if (sequenceQuality != NULL) {
            
            this->setInSequenceQuality(sequenceQuality);
            
            threadLog->add("Segment sequence quality set");
            
        }
        
        this->setACGT(A, C, G, T, N);
        
        threadLog->add("Increased ACGT counts");
        
        this->setLowerCount(lowerCount);

        threadLog->add("Increased total count of lower bases");
        
        seqSize = *A + *C + *G + *T;
        
    }else{
        
        seqSize = *lowerCount;
        
        this->setLowerCount(&seqSize);
        
        threadLog->add("No seq input. Length (" + std::to_string(seqSize) + ") recorded in lower count");
        
    }
    
}
    
void InSegment::setSeqHeader(std::string h) {
    seqHeader = h;
}

void InSegment::setSeqComment(std::string c) {
    seqComment = c;
}

void InSegment::setInSequence(std::string* s) {
    inSequence = s;
}

void InSegment::setInSequenceQuality(std::string* q) {
    inSequenceQuality = q;
}

void InSegment::setSeqTags(std::vector<Tag>* t) {
    tags = *t;
}

void InSegment::setuId(unsigned int i) { // absolute id
    uId = i;
}

void InSegment::setiId(unsigned int i) { // temporary id, internal to scaffold
    iId = i;
}

void InSegment::setSeqPos(unsigned int i) { // temporary id, internal to scaffold
    seqPos = i;
}

std::string InSegment::getSeqHeader() {
    return seqHeader;
}

std::string InSegment::getSeqComment() {
    return seqComment;
}

std::vector<Tag> InSegment::getTags() {
    return tags;
}

std::string InSegment::getInSequence(unsigned int start, unsigned int end) {
    
    if (inSequence == NULL) {
        
        return "*";
        
    }else{
    
        return start != 0 || end != 0 ? inSequence->substr(start-1, end-start+1) : *inSequence;
        
    }
    
}

std::string* InSegment::getInSequencePtr() {
   
    return inSequence;
    
}

std::string InSegment::getInSequenceQuality(unsigned int start, unsigned int end) {
    
    if (inSequenceQuality != NULL) {
    
        return start != 0 || end != 0 ? inSequenceQuality->substr(start-1, end-start+1) : *inSequenceQuality;
        
    }else{
        
        return "";
        
    }
    
}

unsigned int InSegment::getSeqPos() {
    
    return seqPos;

}

uint64_t InSegment::getSegmentLen(uint64_t start, uint64_t end) {
    
    if (inSequence == NULL) {
        
        return lowerCount;
        
    }else{
    
        return start != 0 || end != 0 ? end-start+1 : inSequence->size(); // need to sum long long int to prevent size() overflow
        
    }
    
}

unsigned int InSegment::getuId() { // absolute id
    
    return uId;
}

unsigned int InSegment::getsId() { // temporary id, internal to scaffold
    
    return iId;
}

void InSegment::setACGT(uint64_t* a, uint64_t* c, uint64_t* g, uint64_t* t, uint64_t* n) {
    
    A = *a;
    C = *c;
    G = *g;
    T = *t;
    if (n != NULL)
        N = *n;
    
}

void InSegment::setLowerCount(uint64_t* C) {
    
    lowerCount = *C;
    
}

uint64_t InSegment::getA() {
    
    return A;
}

uint64_t InSegment::getC() {
    
    return C;
}

uint64_t InSegment::getG() {
    
    return G;
}

uint64_t InSegment::getT() {
    
    return T;
}


uint64_t InSegment::getN() {

    return N;
}

unsigned int InSegment::getLowerCount(uint64_t start, uint64_t end) {
    
    if (start == 0 || end == 0) {
        
        return lowerCount;
        
    }else{
        
        uint64_t lowerCountSubset = 0;
        
        for (char base : *inSequence) { // need to fix this loop
            
            if (islower(base)) {
                
                ++lowerCountSubset;
                
            }
            
        }
        
        return lowerCountSubset;
        
    }

}

double InSegment::computeGCcontent() {
    
    double GCcontent = (double) (G + C) / (G + C + A + T) * 100;
    
    return GCcontent;
}

bool InSegment::trimSegment(uint64_t start, uint64_t end) {
    
    if (start > end) {
        std::cerr<<"Trim segment start ("<<+start<<") > end ("<<+end<<"). Exiting.\n";
        exit(EXIT_FAILURE);
    }
    
    inSequence->erase(start, end-start);
    
    if (inSequenceQuality != NULL)
        inSequenceQuality->erase(start, end-start);
    
    return true;
}

bool InSegment::updateSegmentCounts(uint64_t start, uint64_t end) {
    
        std::string newSeq = inSequence->substr(start, end-start);
    
        for(char& base : newSeq) {
    
            switch (base) {
                case 'A':
                case 'a':{
    
                    A--;
                    break;
    
                }
                case 'C':
                case 'c':{
    
                    C--;
                    break;
    
                }
                case 'G':
                case 'g': {
    
                    G--;
                    break;
    
                }
                case 'T':
                case 't': {
    
                    T--;
                    break;
    
                }
    
            }
    
        }
    
    return true;
    
}

bool InSegment::rvcpSegment() {

    *inSequence = revCom(*inSequence);

    return true;
    
}

bool InSegment::invertSegment() {

    *inSequence = rev(*inSequence);
    
    if (inSequenceQuality != NULL) {
    
        *inSequenceQuality = rev(*inSequenceQuality);
    
    }
        
    return true;
    
}

bool InSegment::isCircular(std::vector<unsigned int>* circularSegments) {
    
    return std::binary_search(circularSegments->begin(), circularSegments->end(), this->uId);
    
}

char* InSegment::first() {
    
    return &inSequence->front();
    
}

// GAPS
void InGap::newGap(unsigned int uId, unsigned int sId1, unsigned int sId2, const char& sId1or, const char& sId2or, uint64_t dist, std::string gHeader, std::vector<Tag> tags) {
    
    this->gHeader = gHeader;
    this->uId = uId;
    this->sId1 = sId1;
    this->sId2 = sId2;
    this->sId1Or = sId1or;
    this->sId2Or = sId2or;
    this->dist = dist;
    this->tags = tags;
    
}

void InGap::setuId(unsigned int i) { // absolute id
    uId = i;
}

void InGap::setiId(unsigned int i) { // temporary id, internal to scaffold
    iId = i;
}

void InGap::setsId1(unsigned int i) {
    sId1 = i;
}

void InGap::setsId2(unsigned int i) {
    sId2 = i;
}

void InGap::setDist(uint64_t i) {
    dist = i;
}

std::string InGap::getgHeader() {
    
    return gHeader;
    
}

unsigned int InGap::getuId() {
    
    return uId;
    
}

unsigned int InGap::getsId1() {
    
    return sId1;
    
}

char InGap::getsId1Or() {
    
    return sId1Or;
    
}

unsigned int InGap::getsId2() {
    
    return sId2;
    
}

char InGap::getsId2Or() {
    
    return sId2Or;
    
}

uint64_t InGap::getDist(uint64_t start, uint64_t end) {
    
    return start != 0 || end != 0 ? end-start+1 : dist;
    
}

std::vector<Tag> InGap::getTags() {
    
    return tags;
    
}

void InEdge::newEdge(unsigned int eUId, unsigned int sId1, unsigned int sId2, const char& sId1Or, const char& sId2Or, std::string cigar, std::string eHeader, std::vector<Tag> tags) {
    
    this->eUId = eUId;
    this->sId1 = sId1;
    this->sId2 = sId2;
    this->sId1Or = sId1Or;
    this->sId2Or = sId2Or;
    this->cigar = cigar;
    this->eHeader = eHeader;
    this->tags = tags;
    
}

bool InEdge::operator==(const InEdge& e) const {
    return sId1 == e.sId1 && sId2 == e.sId2 && sId1Or == e.sId1Or && e.sId2Or == sId2Or;
}

void InEdge::seteUId(unsigned int i) { // absolute id
    eUId = i;
}

void InEdge::seteId(unsigned int i) { // temporary id, internal to scaffold
    eId = i;
}

void InEdge::setsId1(unsigned int i) {
    sId1 = i;
}

void InEdge::setsId2(unsigned int i) {
    sId2 = i;
}

void InEdge::setSeqTags(std::vector<Tag>* t) {
    tags = *t;
}

void InEdge::appendTag(Tag t) {
    
    tags.push_back(t);
    
}

std::string InEdge::geteHeader() {
    
    return eHeader;
    
}

std::string InEdge::getCigar() {
    
    return cigar;
    
}

unsigned int InEdge::geteUId() {
    
    return eUId;
    
}

unsigned int InEdge::geteId() {
    
    return eId;
    
}

unsigned int InEdge::getsId1() {
    
    return sId1;
    
}

char InEdge::getsId1Or() {
    
    return sId1Or;
    
}

unsigned int InEdge::getsId2() {
    
    return sId2;
    
}

char InEdge::getsId2Or() {
    
    return sId2Or;

}

std::vector<Tag> InEdge::getTags() {
    return tags;
}

bool InEdge::isSelf() {
    return this->sId1 == this->sId2 ? true : false;
}
    
void InPath::newPath(unsigned int pUid, std::string h, std::string c, unsigned int seqpos) {
    
    pHeader = h;
    pComment = c;
    pathComponents.clear();
    pUId = pUid;
    seqPos = seqpos;

}

void InPath::setpUId(unsigned int pUid) {
    
    pUId = pUid;

}

void InPath::setHeader(std::string pheader) {
    
    pHeader = pheader;

}

void InPath::setComment(std::string c) {
    pComment = c;
}

void InPath::add(PathType type, unsigned int UId, char sign, uint64_t start, uint64_t end) {
    
    pathComponents.push_back({type, UId, sign, start, end});
    
}

void InPath::append(std::vector<PathComponent> components) {

    pathComponents.insert(std::end(pathComponents), std::begin(components), std::end(components));
    
}

void InPath::insert(std::vector<PathComponent>::iterator prevComponent, PathComponent newComponent) {

    pathComponents.insert(prevComponent, newComponent);
    
}

void InPath::clearPath() {
    
    pathComponents.clear();
    
}

void InPath::setComponents(std::vector<PathComponent> newComponents) {

    pathComponents = newComponents;
    
}

std::vector<PathComponent> InPath::getComponents() {
    
    return pathComponents;
    
}

std::vector<PathComponent>* InPath::getComponentsByRef() {
    
    return &pathComponents;
    
}

unsigned int InPath::getpUId() {
    
    return pUId;
    
}

std::string InPath::getHeader() {
    
    return pHeader;
    
}

std::string InPath::getComment() {
    
    return pComment;

}

unsigned int InPath::getSeqPos() {
    
    return seqPos;

}

unsigned int InPath::getContigN() {
    
    return contigN;
    
}

uint64_t InPath::getLen() {
    
    return length;
    
}

uint64_t InPath::getA() {
    
    return A;
    
}

uint64_t InPath::getC() {
    
    return C;
    
}

uint64_t InPath::getG() {
    
    return G;
    
}

uint64_t InPath::getT() {
    
    return T;
    
}

uint64_t InPath::getSegmentLen() {
    
    return segmentLength;
    
}

uint64_t InPath::getLowerCount() {
    
    return lowerCount;
    
}

void InPath::revCom() {
    
    revComPathComponents(pathComponents);

}

void InPath::increaseContigN() {
    
    contigN++;

}

void InPath::increaseGapN() {
    
    contigN++;

}

void InPath::increaseLen(uint64_t n) {
    
    length += n;

}

void InPath::increaseSegmentLen(uint64_t n) {
    
    segmentLength += n;

}

void InPath::increaseLowerCount(uint64_t n) {
    
    lowerCount += n;

}

void InPath::increaseA(uint64_t n) {
    
    A += n;

}

void InPath::increaseC(uint64_t n) {
    
    C += n;

}

void InPath::increaseG(uint64_t n) {
    
    G += n;

}

void InPath::increaseT(uint64_t n) {
    
    T += n;

}

void InPath::reinitializeCounts() {
    
    length = 0, segmentLength = 0, lowerCount = 0, A = 0, C = 0, G = 0, T = 0;

}

void InSegment::addVariants(std::vector<std::vector<DBGpath>> variants) {
    this->variants = variants;
}

std::vector<std::vector<DBGpath>>& InSegment::getVariants() {
    return variants;
}
