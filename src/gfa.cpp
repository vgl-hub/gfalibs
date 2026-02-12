#include <stdlib.h>
#include <string>
#include <vector>

#include <parallel-hashmap/phmap.h>

#include "log.h"
#include "global.h"
#include "bed.h"
#include "struct.h"
#include "functions.h" // global functions
#include "gfa-lines.h"
#include "threadpool.h"
#include "uid-generator.h"

#include "gfa.h"

InSequences::~InSequences() {
    for (InSegment* p : inSegments)
        delete p;    
}

InGap InSequences::pushbackGap(Log* threadLog, InPath* path, std::string* seqHeader, unsigned int* iId, uint64_t &dist, char sign, unsigned int uId1, unsigned int uId2) {
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    threadLog->add("Processing gap " + *seqHeader+"." + std::to_string(*iId) + " (uId: " + std::to_string(uId.get()) + ", iId: " + std::to_string(*iId) + ")");
    
    InGap gap;
    
    gap.newGap(uId.get(), uId1, uId2, '+', sign, dist, *seqHeader+"."+std::to_string(*iId));
    
    insertHash(*seqHeader+"."+std::to_string(*iId), uId.get());
    
    path->add(GAP, uId.get(), '0');
    
    (*iId)++; // number of gaps in the current scaffold
    uId.next(); // unique numeric identifier
    
    lck.unlock();
    
    dist=0;
    
    return gap;
    
}

InSegment* InSequences::pushbackSegment(unsigned int currId, Log* threadLog, InPath* path, std::string* seqHeader, std::string* seqComment, std::string* sequence, unsigned int* iId, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount, uint32_t seqPos, uint64_t sStart, uint64_t sEnd, std::string* sequenceQuality) {
    
    std::string* sequenceSubSeq = new std::string;
    
    *sequenceSubSeq = sequence->substr(sStart, sEnd + 1 - sStart);
    
    if (sequenceQuality != NULL) {
        
        std::string* sequenceQualitySubSeq = new std::string;
        
        *sequenceQualitySubSeq = sequenceQuality->substr(sStart, sEnd + 1 - sStart);
        
        sequenceQuality = sequenceQualitySubSeq;
        
    }
    
    InSegment* inSegment = new InSegment;
    
    inSegment->set(threadLog, currId, *iId, *seqHeader+"."+std::to_string(*iId), seqComment, sequenceSubSeq, A, C, G, T, lowerCount, seqPos, sequenceQuality);
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    insertHash(*seqHeader+"."+std::to_string(*iId), currId);
    
    lck.unlock();
    
    path->add(SEGMENT, currId, '+');
    
    (*iId)++; // number of segments in the current scaffold
    
    *A = 0, *C = 0, *G = 0, *T = 0, *lowerCount = 0;
    
    return inSegment;
    
}

bool InSequences::traverseInSequence(Sequence* sequence, int hc_cutoff) { // traverse the sequence to split at gaps and measure sequence properties
    
    Log threadLog;
    threadLog.setId(sequence->seqPos);
    std::vector<std::pair<uint64_t, uint64_t>> bedCoords;
    if(hc_cutoff != -1)
        homopolymerCompress(sequence->sequence, bedCoords, hc_cutoff);
    
    std::vector<InSegment*> newSegments;
    std::vector<InGap> newGaps;
    uint64_t pos = 0, // current position in sequence
    hc_index=0, // used with homopolymer compression
    A = 0, C = 0, G = 0, T = 0,
    dist = 0, // gap size
    lowerCount = 0,
	sStart = 0, sEnd = 0, // segment start and end
	count = 1; // hc;
	
    unsigned int
    currId = 0, nextId = 0, // temporarily store the id of a segment to connect gaps
	iId = 1; // scaffold feature internal identifier
    char sign = '+';
    bool wasN = false;
    
    sequence->sequence->erase(std::remove(sequence->sequence->begin(), sequence->sequence->end(), '\n'), sequence->sequence->end());
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    lck.lock();
    InPath path;
    
    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find (sequence->header); // get the headers to uIds table to look for the header

    if (got == headersToIds.end()) // this is the first time we see this path name
        insertHash(sequence->header, uId.get());
    else {
        fprintf(stderr, "Error: path name already exists (%s). Terminating.\n", sequence->header.c_str());
        exit(1);
    }
    
    path.newPath(uId.get(), sequence->header, "", sequence->seqPos);
    threadLog.add("Processed sequence: " + sequence->header + " (uId: " + std::to_string(uId.get()) + ")");
    uId.next();
    currId = uId.get();
    uId.next();
    lck.unlock();
    
    if (sequence->comment != "")
        path.setComment(sequence->comment);

    uint64_t seqLen = sequence->sequence->size()-1;
    
    for (char &base : *sequence->sequence) {

        count = 1; // GF, this functionality added by AB is conceptually wrong and should be fixed: it will return the original counts for ACGTN but everything else won't be affected (eg total length will be in hom-compressed space). The hashtable could be further used to output all homopolymer locations
        if(hc_cutoff != -1 && hc_index < bedCoords.size() && pos == bedCoords[hc_index].first) {
            count = bedCoords[hc_index].second - bedCoords[hc_index].first;
            ++hc_index;
        }

        if (islower(base))
            lowerCount+=count;

        switch (base) {

            case 'N':
            case 'n':
            case 'X':
            case 'x': {

                dist+=count;

                if (!wasN && pos>0) { // gap start and gap not at the start of the sequence

                    sEnd = pos - 1;
                    newSegments.push_back(pushbackSegment(currId, &threadLog, &path, &sequence->header, &sequence->comment, sequence->sequence, &iId, &A, &C, &G, &T, &lowerCount, sequence->seqPos, sStart, sEnd, sequence->sequenceQuality));
                    
                    lck.lock();
                    uId.next();
                    lck.unlock();
                }
                if(pos == seqLen) { // end of scaffold, terminal gap

                    sign = '-';
                    newGaps.push_back(pushbackGap(&threadLog, &path, &sequence->header, &iId, dist, sign, currId, currId));
                }
                wasN = true;
                break;
            }
            default: {

                switch (base) {
                    case 'A':
                    case 'a':{
                        A+=count;
                        break;
                    }
                    case 'C':
                    case 'c':{
                        C+=count;
                        break;
                    }
                    case 'G':
                    case 'g': {
                        G+=count;
                        break;
                    }
                    case 'T':
                    case 't': {
                        T+=count;
                        break;
                    }
                }
                if (wasN) { // internal gap end
                    
                    lck.lock();
                    uId.next();
                    nextId = uId.get();
                    uId.next();
                    lck.unlock();
                    
                    if (newSegments.size() == 0)
                        currId = nextId;

                    sStart = pos;
                    newGaps.push_back(pushbackGap(&threadLog, &path, &sequence->header, &iId, dist, sign, currId, nextId));
                    currId = nextId;
                }

                if (pos == seqLen) {

                    sEnd = pos;
                    newSegments.push_back(pushbackSegment(currId, &threadLog, &path, &sequence->header, &sequence->comment, sequence->sequence, &iId, &A, &C, &G, &T, &lowerCount, sequence->seqPos, sStart, sEnd, sequence->sequenceQuality));
                    lck.lock();
                    uId.next();
                    lck.unlock();
                }
                wasN = false;
            }
        }
        pos++;
    }
    delete sequence;
    lck.lock();
    inGaps.insert(std::end(inGaps), std::begin(newGaps), std::end(newGaps));
    threadLog.add("Segments added to segment vector");
    inSegments.insert(std::end(inSegments), std::begin(newSegments), std::end(newSegments));
    threadLog.add("Gaps added to segment vector");
    inPaths.push_back(path);
    threadLog.add("Added fasta sequence as path");
    logs.push_back(threadLog);
    lck.unlock();
    return true;
}

bool InSequences::traverseInSegmentWrapper(Sequence* sequence, std::vector<Tag> inSequenceTags) {
    
    traverseInSegment(sequence, inSequenceTags);
    return true;
    
}

InSegment* InSequences::traverseInSegment(Sequence* sequence, std::vector<Tag> inSequenceTags, uint32_t sId) { // traverse the segment
    
    Log threadLog;
    
    threadLog.setId(sequence->seqPos);
    
    uint64_t A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
    unsigned int sUId = 0;
    
    for (char &base : *sequence->sequence) {
        
        if (islower(base)) {
            
            lowerCount++;
            
        }
                
        switch (base) {
            case 'A':
            case 'a':{
                
                A++;
                break;
                
            }
            case 'C':
            case 'c':{
                
                C++;
                break;
                
            }
            case 'G':
            case 'g': {
                
                G++;
                break;
                
            }
            case 'T':
            case 't': {
                
                T++;
                break;
                
            }
                
            case '*': {
                
                auto tag = find_if(inSequenceTags.begin(), inSequenceTags.end(), [](Tag& obj) {return checkTag(obj.label, "LN");}); // find if length tag is present in the case sequence is missing
                
                if (tag != inSequenceTags.end()) {
                    
                    lowerCount = stol(tag->content);
                    
                }
                    
                break;
                
            }
                
        }
            
    }
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find(sequence->header); // get the headers to uIds table to look for the header
    
    if (got == headersToIds.end()) { // this is the first time we see this segment
        sUId = uId.next();
        insertHash(sequence->header, sUId);
    }else{
        sUId = got->second;
    }
    
    InSegment* inSegment = new InSegment;
    
    inSegment->set(&threadLog, sUId, sId, sequence->header, &sequence->comment, sequence->sequence, &A, &C, &G, &T, &lowerCount, sequence->seqPos, sequence->sequenceQuality, &inSequenceTags);
            
    inSegments.push_back(inSegment);
    
    logs.push_back(threadLog);
    
    lck.unlock();
    
    return inSegment;
    
}

void InSequences::appendSequence(Sequence* sequence, int hc_cutoff) { // method to append a new sequence from a fasta
        
    threadPool.queueJob([=]{ return traverseInSequence(sequence, hc_cutoff); });
    
    if(verbose_flag) {std::cerr<<"\n";};
    
}

void InSequences::appendSegment(Sequence* sequence, std::vector<Tag> inSequenceTags) { // method to append a new segment from a gfa
    
    lg.verbose("Segment read");
    
    threadPool.queueJob([=]{ return traverseInSegmentWrapper(sequence, inSequenceTags); });
    
    if(verbose_flag) {std::cerr<<"\n";};
    
}

InSegment* InSequences::getInSegment(unsigned int sId) {
    
    auto inSegment = find_if(inSegments.begin(), inSegments.end(), [sId](InSegment* obj) {return obj->getuId() == sId;}); // given a uId, find it in nodes
    
    if (inSegment == inSegments.end()) {
    
        fprintf(stderr, "Error: could not find segment (sId: %i).\n", sId); exit(1);
        
    }
        
    return *inSegment;
    
}

std::vector<InSegment*>* InSequences::getInSegments() {
    
    return &inSegments;
    
}

std::vector<InGap>* InSequences::getInGaps() {
    
    return &inGaps;
    
}

InPath InSequences::getInPath(unsigned int pId) {
    
    auto inPath = find_if(inPaths.begin(), inPaths.end(), [pId](InPath& obj) {return obj.getpUId() == pId;}); // given a uId, find it in nodes
    
    if (inPath == inPaths.end()) {
    
        fprintf(stderr, "Error: could not find path (pId: %i).\n", pId); exit(1);
        
    }
        
    return *inPath;
    
}

std::vector<InPath> InSequences::getInPaths() {
    
    return inPaths;
    
}

uint64_t InSequences::getTotScaffLen() {
    
    return totScaffLen;
    
}

uint64_t InSequences::getTotSegmentLen() {
    
    for (std::vector<InSegment*>::iterator inSegment = inSegments.begin(); inSegment != inSegments.end(); inSegment++) {
        
        totSegmentLen += (*inSegment)->getA() + (*inSegment)->getC() + (*inSegment)->getG() + (*inSegment)->getT();
        
    }

    return totSegmentLen;
    
}

void InSequences::changeTotGapLen(int64_t gapLen) {
    totGapLen += gapLen;
}

unsigned int InSequences::getGapNScaffold() {
    
    return gapLens.size();
    
}

unsigned int InSequences::getGapN() {
    
    return inGaps.size();
    
}

unsigned int InSequences::getEdgeN() {
    
    return inEdges.size();
    
}

unsigned int InSequences::getPathN() {
    
    return inPaths.size();
    
}

void InSequences::recordScaffLen(uint64_t seqLen) {
    
    scaffLens.push_back(seqLen);
    
}

void InSequences::recordGapLen(uint64_t gapLen) {
    gapLens.push_back(gapLen);
}

void InSequences::evalNstars(char type, uint64_t gSize) { // switch between scaffold, contig, gap while computing N* statistics
    
    switch(type) {
            
        case 's': {
            
            computeNstars(scaffLens, scaffNstars, scaffLstars, &scaffNGstars, &scaffLGstars, gSize);
            break;
            
        }
            
        case 'c': {
            
            computeNstars(contigLens, contigNstars, contigLstars, &contigNGstars, &contigLGstars, gSize);
            break;
            
        }
            
        case 'g': {
            
            computeNstars(gapLens, gapNstars, gapLstars);
            break;
            
        }
            
    }
    
}

void InSequences::evalAuN(char type, uint64_t gSize) { // switch between scaffold, contig, gap while computing N* statistics
    
    switch(type) {
            
        case 's': {
            
            computeAuN(scaffLens, scaffAuN, &scaffAuNG, gSize);
            break;
            
        }
            
        case 'c': {
            
            computeAuN(contigLens, contigAuN, &contigAuNG, gSize);
            break;
            
        }
            
        case 'g': {
            
            computeAuN(gapLens, gapAuN);
            break;
            
        }
            
    }
    
}

void InSequences::computeAuN(std::vector<uint64_t>& lens, double& auN, double* auNG, uint64_t gSize) {// compute N* statistics
    
    uint64_t totLen = 0;
    
    for(std::vector<uint64_t>::iterator it = lens.begin(); it != lens.end(); ++it) // find total length
        totLen += *it;
    
    for(unsigned int i = 0; i < lens.size(); i++) { // for each length
        
        auN += (double) lens[i] * lens[i] / totLen; // the area under the curve is the length (height) * fraction of the total length
        
        if (gSize > 0) {
            
            *auNG += (double) lens[i] * lens[i] / gSize;
            
        }
        
    }
    
}

unsigned int InSequences::getScaffN() {
    
    return scaffN;
    
}

unsigned int InSequences::getTotContigN() {
    
    return contigN;
    
}

std::vector <uint64_t> InSequences::getScaffNstars() {
    
    return scaffNstars;
    
}

std::vector <uint64_t> InSequences::getScaffNGstars() {
    
    return scaffNGstars;
    
}

std::vector <unsigned int> InSequences::getScaffLstars() {
    
    return scaffLstars;
    
}

std::vector <unsigned int> InSequences::getScaffLGstars() {
    
    return scaffLGstars;
    
}

std::vector <uint64_t> InSequences::getContigNstars() {
    
    return contigNstars;
    
}

std::vector <uint64_t> InSequences::getContigNGstars() {
    
    return contigNGstars;
    
}

std::vector <unsigned int> InSequences::getContigLstars() {
    
    return contigLstars;
    
}

std::vector <unsigned int> InSequences::getContigLGstars() {
    
    return contigLGstars;
    
}

std::vector <uint64_t> InSequences::getGapNstars() {
    
    return gapNstars;
    
}

std::vector <unsigned int> InSequences::getGapLstars() {
    
    return gapLstars;
    
}

uint64_t InSequences::getScaffN50() {
    
    return scaffNstars[4];
    
}

uint64_t InSequences::getScaffNG50() {
    
    return scaffNGstars[4];
    
}

unsigned int InSequences::getScaffL50() {
    
    return scaffLstars[4];
    
}

unsigned int InSequences::getScaffLG50() {
    
    return scaffLGstars[4];
    
}

double InSequences::getScaffauN() {
    
    return scaffAuN;
    
}

double InSequences::getScaffauNG() {
    
    return scaffAuNG;
    
}

double InSequences::getContigauN() {
    
    return contigAuN;
    
}

double InSequences::getContigauNG() {
    
    return contigAuNG;
    
}

double InSequences::getGapauN() {
    
    return gapAuN;
    
}

unsigned int InSequences::getSegmentN() {
    
    return inSegments.size();
    
}

uint64_t InSequences::getContigN50() {
    
    return contigNstars[4]; // middle value
    
}

unsigned int InSequences::getContigNG50() {
    
    return contigNGstars[4];
    
}

unsigned int InSequences::getContigL50() {
    
    return contigLstars[4];
    
}

unsigned int InSequences::getContigLG50() {
    
    return contigLGstars[4];
    
}

unsigned int InSequences::getGapN50() {
    
    return gapNstars[4];
    
}

unsigned int InSequences::getGapL50() {
    
    return gapLstars[4];
    
}

uint64_t InSequences::getLargestScaffold() {
    
    return scaffLens.size() == 0 ? 0 : scaffLens[0]; // sorted during N/L* computation
    
}

uint64_t InSequences::getSmallestScaffold() {
    
    return scaffLens.size() == 0 ? 0 : scaffLens.back(); // sorted during N/L* computation
    
}


uint64_t InSequences::getLargestContig() {
    
    return contigLens.size() == 0 ? 0 : contigLens[0]; // sorted during N/L* computation
    
}

uint64_t InSequences::getSmallestContig() {
    
    return contigLens.size() == 0 ? 0 : contigLens.back(); // sorted during N/L* computation
    
}


unsigned int InSequences::getLargestGap() {
    
    return gapLens.size() == 0 ? 0 : gapLens[0]; // sorted during N/L* computation
    
}

unsigned int InSequences::getSmallestGap() {
    
    return gapLens.size() == 0 ? 0 : gapLens.back(); // sorted during N/L* computation
    
}

double InSequences::computeAvgScaffLen() {
    
    return (double) InSequences::totScaffLen/scaffN;
    
}

uint64_t InSequences::getTotContigLen () {
    
    totContigLen = 0;
    
    for (std::vector<uint64_t>::iterator contigLen = contigLens.begin(); contigLen != contigLens.end(); contigLen++)
        totContigLen += *contigLen;

    return totContigLen;
}

double InSequences::computeAvgContigLen() {
    
    return (double) getTotContigLen()/contigLens.size();
    
}

double InSequences::computeAvgSegmentLen() {
    
    return (double) totSegmentLen/inSegments.size();
    
}

uint64_t InSequences::getTotGapLen() {
    
    totGapLen = 0;
    
    for (std::vector<uint64_t>::iterator gapLen = gapLens.begin(); gapLen != gapLens.end(); gapLen++) {
        
        totGapLen += *gapLen;
        
    }
    
    return totGapLen;
    
}

double InSequences::computeAverageGapLen() {
    
    return totGapLen == 0 ? 0 : (double) getTotGapLen()/gapLens.size();
    
}

uint64_t InSequences::getTotA() {
    
    return totA;
}

uint64_t InSequences::getTotC() {
    
    return totC;
}

uint64_t InSequences::getTotG() {
    
    return totG;
}

uint64_t InSequences::getTotT() {
    
    return totT;
}

uint64_t InSequences::getTotLowerCount() {
    
    return totLowerCount;
}

double InSequences::computeGCcontent() {
    
    double GCcontent;
    
    if (inSegments.size()>0) {
        
        GCcontent = (double) (totC + totG) / (totA + totC + totG + totT) * 100;
        
    }else{
        
        GCcontent = 0;
        
    }
    
    return GCcontent;
}

//gfa methods
bool InSequences::addGap(InGap inGap) {
    
    inGaps.push_back(inGap);

    lg.verbose("Gap added to gap vector");
    
    return true;
    
}

bool InSequences::addPath(InPath path) {
    
    inPaths.push_back(path);

    lg.verbose("Path added to path vector");
    
    pathN++;
    
    lg.verbose("Increased path counter");
    
    return true;
    
}

std::vector<InGap> InSequences::getGaps() { // return gfa gaps vector
    
    return inGaps;
    
}

std::vector<InEdge>* InSequences::getEdges() { // return gfa edge vector
    
    return &inEdges;
    
}

bool InSequences::appendEdge(InEdge edge) {
    
    std::lock_guard<std::mutex> lck(mtx);
    
    if (edge.geteUId() == 0)
        edge.seteUId(uId.next());
    
    inEdges.push_back(edge);
    lg.verbose("Edge added to edge vector");
    return true;
}

//sorting methods

void InSequences::sortSegmentsByOriginal(){
    
    sort(inSegments.begin(), inSegments.end(), [](InSegment* one, InSegment* two){
        
        if(one->getSeqPos() != two->getSeqPos())
          return (one->getSeqPos() < two->getSeqPos());
        return one->getsId() < two->getsId();
        
    });
    
}

void InSequences::sortEdgesByOriginal(){
    
    sort(inEdges.begin(), inEdges.end(), [](InEdge one, InEdge two){
        
        if(one.geteUId() != two.geteUId())
          return (one.geteUId() < two.geteUId());
        return one.geteId() < two.geteId();
        
    });
    
}

void InSequences::sortPathsByOriginal(){
    
    sort(inPaths.begin(), inPaths.end(), [](InPath& one, InPath& two){return one.getSeqPos() < two.getSeqPos();});
    
}

void InSequences::sortPathsByNameAscending(){
    
    sort(inPaths.begin(), inPaths.end(), [](InPath& one, InPath& two){return one.getHeader() < two.getHeader();});
    
}

void InSequences::sortPathsByNameDescending(){
    
    sort(inPaths.begin(), inPaths.end(), [](InPath& one, InPath& two){return one.getHeader() > two.getHeader();});
    
}

void InSequences::sortPathsByList(std::vector<std::string> headerList){
    
    int index1 = 0, index2 = 0;
    auto comp = [&](InPath& one, InPath& two)-> bool { // lambda function for custom sorting
    auto it = find(headerList.begin(), headerList.end(), one.getHeader());
    if (it != headerList.end()) { // if element one was found
        index1 = it - headerList.begin(); // calculating the index
    }else {
        std::cout<<"Error: sequence missing from sort list (" << one.getHeader() << ")\n";
        exit(1);
    }
    it = find(headerList.begin(), headerList.end(), two.getHeader());
    if (it != headerList.end()) { // if element two was found
        index2 = it - headerList.begin(); // calculating the index
    }else {
        std::cout<<"Error: sequence missing from sort list ("<<two.getHeader()<<")\n";
        exit(1);
    }
        return index1<index2;
    };
    sort(inPaths.begin(), inPaths.end(), comp);
}

void InSequences::sortPathsBySize(bool largest){
    
    auto comp = [&](InPath& one, InPath& two)-> bool { // lambda function for custom sorting
    
        std::vector<PathComponent> pathComponents;
        
		unsigned int uId = 0, sIdx = 0, gIdx = 0;
		uint64_t size1 = 0, size2 = 0;
            
        pathComponents = one.getComponents();
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            uId = component->id;
            
            if (component->componentType == SEGMENT) {
            
                auto sId = find_if(inSegments.begin(), inSegments.end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                
                if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index
                
                size1 += inSegments[sIdx]->getInSequence().size();
                
            }else{
                
                auto gId = find_if(inGaps.begin(), inGaps.end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                
                if (gId != inGaps.end()) {gIdx = std::distance(inGaps.begin(), gId);} // gives us the segment index
                
                size1 += inGaps[gIdx].getDist();
                
            }
            
        }
            
        pathComponents = two.getComponents();
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            uId = component->id;
            
            if (component->componentType == SEGMENT) {
                auto sId = find_if(inSegments.begin(), inSegments.end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index
                size2 += inSegments[sIdx]->getInSequence().size();
            }else{
                auto gId = find_if(inGaps.begin(), inGaps.end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                if (gId != inGaps.end()) {gIdx = std::distance(inGaps.begin(), gId);} // gives us the segment index
                size2 += inGaps[gIdx].getDist();
            }
        }
        if(largest)
            return size1<size2;
        else
            return size1>size2;
    };
        
    sort(inPaths.begin(), inPaths.end(), comp);
    
}

//gfa methods
void InSequences::insertHash(const std::string &segHeader, unsigned int i) {
    headersToIds.insert({segHeader, i});
    idsToHeaders.insert({i, segHeader});
}

void InSequences::updateHash(const std::string &segHeader, unsigned int i) {
    headersToIds[segHeader] = i;
    idsToHeaders[i] = segHeader;
}

void InSequences::eraseHash(const std::string &segHeader, unsigned int i) {
    headersToIds.erase(segHeader);
    idsToHeaders.erase(i);
}

unsigned int InSequences::getuId() {
    return uId.get();
}

phmap::flat_hash_map<std::string, unsigned int>* InSequences::getHash1() {
    return &headersToIds;
}

phmap::flat_hash_map<unsigned int, std::string>* InSequences::getHash2() {
    return &idsToHeaders;
}

void InSequences::buildGraph(std::vector<InGap> const& gaps) { // graph constructor
    
    lg.verbose("Started graph construction");
    
    adjListFW.clear();
    adjListBW.clear();
    
    adjListFW.resize(uId.get()); // resize the adjaciency list to hold all nodes
    adjListBW.resize(uId.get()); // resize the adjaciency list to hold all nodes
    
    for (auto &gap: gaps) // add edges to the graph
    {
        
        lg.verbose("Adding forward gap " + std::to_string(gap.uId) + ": " + idsToHeaders[gap.sId1] + "(" + std::to_string(gap.sId1) + ") " + gap.sId1Or + " " + idsToHeaders[gap.sId2] + "(" + std::to_string(gap.sId2) + ") " + gap.sId2Or + " " + std::to_string(gap.dist));
        
        adjListFW.at(gap.sId1).push_back({gap.sId1Or, gap.sId2, gap.sId2Or, gap.dist, gap.uId}); // insert at gap start gap destination, orientations and weight (gap size)

        lg.verbose("Adding reverse gap " + std::to_string(gap.uId) + ": " + idsToHeaders[gap.sId2] + "(" + std::to_string(gap.sId2) + ") " + edge.sId2Or + " " + idsToHeaders[gap.sId1] + "(" + std::to_string(gap.sId1) + ") " + edge.sId2Or + " " + std::to_string(gap.dist));
        
        adjListBW.at(gap.sId2).push_back({gap.sId2Or, gap.sId1, gap.sId1Or, gap.dist, gap.uId}); // undirected graph
        
    }
    
    lg.verbose("Graph built");
    
    visited.clear();
    
}

void InSequences::buildEdgeGraph() { // graph constructor
    
    lg.verbose("Started edge graph construction");
    adjEdgeList.clear();
    adjEdgeList.resize(uId.get()); // resize the adjaciency list to hold all nodes
    
    for (auto &edge: inEdges) { // add edges to the graph
		
		Edge fwEdge = {edge.sId1Or, edge.sId2, edge.sId2Or};
		
		if (find(adjEdgeList.at(edge.sId1).begin(), adjEdgeList.at(edge.sId1).end(), fwEdge) != adjEdgeList.at(edge.sId1).end()) // add edge only if is not already present
			continue;
        
        lg.verbose("Adding edge: " + idsToHeaders[edge.sId1] + "(" + std::to_string(edge.sId1) + ") " + edge.sId1Or + " " + idsToHeaders[edge.sId2] + "(" + std::to_string(edge.sId2) + ") " + edge.sId2Or);

        adjEdgeList.at(edge.sId1).push_back(fwEdge); // insert at edge start gap destination and orientations
        Edge rvEdge {edge.sId2Or == '+' ? '-' : '+', edge.sId1, edge.sId1Or == '+' ? '-' : '+'};
		adjEdgeList.at(edge.sId2).push_back(rvEdge); // assembly are bidirected by definition
    }
    lg.verbose("Graph built");
    visited.clear();
}

std::vector<std::vector<Edge>>& InSequences::getAdjEdgeList() {
    return adjEdgeList;
}

void InSequences::dfsEdges(unsigned int v, uint64_t* componentLength) { // Depth First Search to explore graph connectivity

   visited[v] = true; // mark the current node as visited
   unsigned int sIdx = 0;
   auto sId = find_if(inSegments.begin(), inSegments.end(), [v](InSegment* obj) {return obj->getuId() == v;}); // given a node Uid, find it
   if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index

   if (adjEdgeList.at(v).size() > 1) { // if the vertex has more than one edge

        *componentLength += inSegments[sIdx]->getSegmentLen();
        char sign = adjEdgeList.at(v).at(0).orientation0;
        unsigned int i = 0;

        for(auto edge : adjEdgeList.at(v)) {
        
            i++;

            if(edge.orientation0 != sign){

                lg.verbose("node: " + idsToHeaders[v] + " --> case a: internal node, multiple edges");
                break;

            }else if (i == adjEdgeList.at(v).size()) {

                lg.verbose("node: " + idsToHeaders[v] + " --> case b: single dead end, multiple edges");
                deadEnds += 1;
            }
        sign = edge.orientation0;
        }
    }else if (adjEdgeList.at(v).size() == 1){ // this is a single dead end
        deadEnds += 1;
        *componentLength += inSegments[sIdx]->getSegmentLen();
        lg.verbose("node: " + idsToHeaders[v] + " --> case c: single dead end, single edge");
    }else if(adjEdgeList.at(v).size() == 0){ // disconnected component (double dead end)
        deadEnds += 2;
        disconnectedComponents++;
        lengthDisconnectedComponents += inSegments[sIdx]->getSegmentLen();
        lg.verbose("node: " + idsToHeaders[v] + " --> case d: disconnected component");
    }
    for (auto i: adjEdgeList[v]) { // recur for all forward vertices adjacent to this vertex
       if (!visited[i.id] && !deleted[i.id])
           dfsEdges(i.id, componentLength); // recurse
    }
}

void InSequences::dfsScaffolds(unsigned int v, uint64_t* scaffSize, uint64_t* A, uint64_t* C, uint64_t* G, uint64_t* T, uint64_t* lowerCount) // Depth First Search to explore graph connectivity
{
    
    visited[v] = true; // mark the current node as visited
    unsigned int idx = 0, a = 0, c = 0, g = 0, t = 0;
    
    bool seqRevCom = false, segRevCom = false;
    
    auto it = find_if(inSegments.begin(), inSegments.end(), [&v](InSegment* obj) {return obj->getuId() == v;}); // given a vertex id, search it in the segment vector
    
    if (it != inSegments.end()) {idx = std::distance(inSegments.begin(), it);} // if found, get its index
    
    if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(adjListFW.at(v).at(0).segmentId == adjListBW.at(v).at(0).segmentId) && !backward) { // if the vertex has exactly one forward and one backward connection and they do not connect to the same vertex (internal node)
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case a: internal node, forward direction");
        
        seqRevCom = (adjListFW.at(v).at(0).orientation0 == '+') ? false : true; // check if sequence should be in forward orientation, if not reverse-complement
        
        if (seqRevCom) {
            
            unsigned int tmpA = *A, tmpC = *C;
            
            *A = *T;
            *C = *G;
            *G = tmpC;
            *T = tmpA;
            
        }
        
        segRevCom = (adjListBW.at(v).at(0).orientation0 == '+') ? false : true; // check if vertex should be in forward orientation, if not reverse-complement
        
        backward = false;
        
    }else if (adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 1){ // this is the final vertex without gaps
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case b: end node, forward direction, no final gap");
        
        segRevCom = (adjListBW.at(v).at(0).orientation0 == '+') ? false : true;
        
        backward = true; // reached the end
        
    }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 2){ // this is the final vertex with terminal gap
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case c: end node, forward direction, final gap");
        
        if(adjListBW.at(v).at(0).segmentId != v) { // make sure you are not using the terminal edge to ascertain direction in case it was edited by sak
        
            segRevCom = (adjListBW.at(v).at(0).orientation0 == '+') ? false : true;
         
            *scaffSize += adjListBW.at(v).at(1).dist;
            
        }else{
        
            segRevCom = (adjListBW.at(v).at(1).orientation0 == '+') ? false : true;
        
            *scaffSize += adjListBW.at(v).at(0).dist;
            
        }
        
        backward = true; // reached the end
        
    }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(adjListFW.at(v).at(0).segmentId == adjListBW.at(v).at(0).segmentId) && backward){ // this is an intermediate vertex, only walking back
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case d: intermediate node, backward direction");
        
        segRevCom = (adjListBW.at(v).at(0).orientation0 == '+') ? false : true;
        
        backward = true;
        
    }else if(adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 0){ // disconnected component
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case e: disconnected component");
        
    }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 0){ // this is the first vertex without gaps
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case f: start node, no gaps");
        
        segRevCom = (adjListFW.at(v).at(0).orientation0 == '+') ? false : true;
        
        backward = false;
        
    }else if (adjListFW.at(v).size() == 2 && adjListBW.at(v).size() == 1){ // this is the first vertex with a terminal gap
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case g: start node, start gap");
        
        segRevCom = (adjListFW.at(v).at(0).orientation0 == '+') ? false : true;
        
        *scaffSize += adjListFW.at(v).at(0).dist;
        
        backward = false;
        
    }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && adjListFW.at(v).at(0).segmentId == adjListBW.at(v).at(0).segmentId) { // if the vertex has exactly one forward and one backward connection and they connect to the same vertex (disconnected component with gap)
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case h: disconnected component with gap");
        
        segRevCom = (adjListFW.at(v).at(0).orientation0 == '+') ? false : true;
        
        *scaffSize += adjListFW.at(v).at(0).dist;
        
        backward = false;
        
    }
    
    *scaffSize += inSegments[idx]->getInSequence().size();
    
    if (!segRevCom) {
        
        a = inSegments[idx]->getA();
        c = inSegments[idx]->getC();
        g = inSegments[idx]->getG();
        t = inSegments[idx]->getT();
        
    }else {
        
        a = inSegments[idx]->getT();
        c = inSegments[idx]->getG();
        g = inSegments[idx]->getC();
        t = inSegments[idx]->getA();
        
    }
    
    *A += a;
    *C += c;
    *G += g;
    *T += t;
    
    *lowerCount += inSegments[idx]->getLowerCount();
    
    for (auto i: adjListFW[v]) { // recur for all forward vertices adjacent to this vertex
        
        if (!visited[i.segmentId] && !deleted[i.segmentId]) {
            
            *scaffSize += i.dist;
            
            dfsScaffolds(i.segmentId, scaffSize, A, C, G, T, lowerCount); // recurse
            
        }
    }
    
    for (auto i: adjListBW[v]) { // recur for all backward vertices adjacent to this vertex
        
        if (!visited[i.segmentId] && !deleted[i.segmentId]) {
            
            *scaffSize += i.dist;
            
            dfsScaffolds(i.segmentId, scaffSize, A, C, G, T, lowerCount); // recurse
            
        }
    }
    
}

std::vector<std::vector<Gap>> InSequences::getAdjListFW() {
    
    return adjListFW;
    
}

std::vector<std::vector<Gap>> InSequences::getAdjListBW() {
    
    return adjListBW;
    
}

bool InSequences::getVisited(unsigned int uId) {
    
    return visited[uId];
    
}

bool InSequences::getDeleted(unsigned int uId) {
    
    return deleted[uId];
    
}

bool InSequences::updateStats() {
    
    scaffLens.clear();
    contigLens.clear();
    gapLens.clear();
    
    totScaffLen = 0, scaffN = 0, contigN = 0, totA = 0, totC = 0, totG = 0, totT = 0, totLowerCount = 0;
    
    for (InPath& inPath : inPaths) { // loop through all paths
        
        walkPath(&inPath);
        
        totScaffLen += inPath.getLen();
     
        lg.verbose("Increased total scaffold length");
        
        scaffN++;
        
        lg.verbose("Increased total scaffold N");
        
        contigN += inPath.getContigN();
        
        lg.verbose("Increased total contig N");
        
        recordScaffLen(inPath.getLen());
        
        lg.verbose("Recorded length of sequence: " + std::to_string(inPath.getLen()));
        
        totA += inPath.getA();
        totC += inPath.getC();
        totG += inPath.getG();
        totT += inPath.getT();
        
        lg.verbose("Increased total ACGT counts");
        
        totLowerCount += inPath.getLowerCount();

        lg.verbose("Increased total count of lower bases");
        
    }
    
    lg.verbose("Updated scaffold statistics");
    
    return true;
    
}

bool InSequences::removeTerminalGaps() {
    
    std::vector<InPath>::iterator pathIt = inPaths.begin(); // first, remove the gaps from the paths they occur in
    std::vector<PathComponent> pathComponents;
    
    unsigned int uId = 0, gIdx = 0;
    
    while (pathIt != end(inPaths)) {
        
        pathComponents = pathIt->getComponents();
        
        for (std::vector<PathComponent>::iterator componentIt = pathComponents.begin(); componentIt != pathComponents.end();) {
            
            if (componentIt->componentType == GAP) {
                
                uId = componentIt->id;
                
                auto gId = find_if(inGaps.begin(), inGaps.end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a gap Uid, find it
                
                if (gId->getsId1() == gId->getsId2()) { // terminal gaps connect on themeselves
                    
                    gIdx = std::distance(pathComponents.begin(), componentIt); // gives us the gap index
                
                    pathIt->pathComponents.erase(pathIt->pathComponents.begin()+gIdx); // remove the gap component from the path
                    
                    lg.verbose("Removed gap from paths");
                    
                    inGaps.erase(gId); // remove the element by position, considering elements that were already removed in the loop

                    changeTotGapLen(-gId->getDist()); // update length of gaps
                    
                    lg.verbose("Removed gap from gap vector");
                    
                }
                
            }
            
            componentIt++;
            
        }
        
        pathIt++;
        
    }
    
    gapLens.clear();
    
    for (unsigned int i = 0; i != inGaps.size(); ++i) { // update statistics to reflect the removal
    
        recordGapLen(inGaps[i].getDist());
        
    }
    
    lg.verbose("Recorded length of gaps in sequence");
    
    return true;
    
}

unsigned int InSequences::getDeadEnds() {

    return deadEnds;

}

unsigned int InSequences::getDisconnectedComponents() {

    return disconnectedComponents;

}

unsigned int InSequences::getLengthDisconnectedComponents() {

    return lengthDisconnectedComponents;

}

// instruction methods

std::vector<InGap> InSequences::getGap(std::string* contig1, std::string* contig2) { // if two contigs are provided, returns all gaps connecting them, if only one contig is provided returns all gaps where it appears

    std::vector<InGap> gaps;
    
    unsigned int sUId1 = headersToIds[*contig1];
    
    if (contig2 != NULL) {
    
        unsigned int sUId2 = headersToIds[*contig2];
    
        auto gId = find_if(inGaps.begin(), inGaps.end(), [sUId1, sUId2](InGap& obj) {return ( // given two vertex ids, search the gap that connects them
            
            (obj.getsId1() == sUId1 && obj.getsId2() == sUId2) || // fw orientation
            (obj.getsId1() == sUId2 && obj.getsId2() == sUId1)    // rv orientation
                                                                                             
        );});
        
        gaps.push_back(*gId);
        
    }else{
        
        for (InGap inGap : inGaps) {
            
            if (inGap.getsId1() == sUId1 || inGap.getsId2() == sUId1) {
                
                gaps.push_back(inGap);
                
            }
            
        }
        
    }
    
    return gaps;
    
}

std::vector<unsigned int> InSequences::removeGaps(std::string* contig1, std::string* contig2) { // if two contigs are provided, remove all gaps connecting them, if only one contig is provided remove all gaps where it appears
    
    std::vector<unsigned int> guIds;

    unsigned int sUId1 = headersToIds[*contig1], gIdx = 0;
    
    if (contig2 != NULL) {
    
        unsigned int sUId2 = headersToIds[*contig2];
    
        auto gId = find_if(inGaps.begin(), inGaps.end(), [sUId1, sUId2](InGap& obj) {return ( // given two vertex ids, search the gap that connects them
            
            (obj.getsId1() == sUId1 && obj.getsId2() == sUId2) || // fw orientation
            (obj.getsId1() == sUId2 && obj.getsId2() == sUId1)    // rv orientation
                                                                                             
        );});
    
        if (gId != inGaps.end()) {
            
            gIdx = std::distance(inGaps.begin(), gId); // gives us the gap index
            
            guIds.push_back((*gId).getuId());
        
            inGaps.erase(inGaps.begin()+gIdx); // remove the element by position, considering elements that were already removed in the loop
            
            changeTotGapLen(-gId->getDist()); // update length of gaps
            
        }else{
            
            fprintf(stderr, "Error: could not find gap between segments (gId1: %s, gId2: %s).\n", contig1->c_str(), contig2->c_str()); exit(1);
            
        }
        
    }else{
        
        auto it = inGaps.begin();
        
        while (it != end(inGaps)) {

            auto gId = find_if(it, inGaps.end(), [sUId1](InGap& obj) {return obj.getsId1() == sUId1 || obj.getsId2() == sUId1;}); // check whether an edge containing the node was found

            if (gId != inGaps.end()) {
                
                gIdx = std::distance(inGaps.begin(), gId); // gives us the gap index
                
                guIds.push_back((*gId).getuId());
        
                inGaps.erase(inGaps.begin()+gIdx); // remove the element by position, considering elements that were already removed in the loop
                
                changeTotGapLen(-gId->getDist()); // update length of gaps
                
            }
            
            it = gId;
            
        }
        
    }
    
    return guIds;
    
}

bool InSequences::flagDeletedSegment(std::string* contig1) { // flag segment as deleted
    
    if (contig1 != NULL) {
        
        unsigned int sIdx = 0, sUId = headersToIds[*contig1];
    
        auto sId = find_if(inSegments.begin(), inSegments.end(), [sUId](InSegment* obj) {return obj->getuId() == sUId;}); // given a node uId, find it
    
        if (sId != inSegments.end()) {sIdx = std::distance(inSegments.begin(), sId);} // gives us the segment index
    
        deleted[sIdx] = true;
        
    }else{
        
        lg.verbose("Cannot detect node: " + *contig1);
        
    }
    
    return true;
    
}

bool InSequences::deleteSegment(std::string sHeader) { // fully delete segment from inSequence object
    
    std::unique_lock<std::mutex> lck (mtx);
    
    uint32_t sIdx, sUId = headersToIds[sHeader];
    auto sId = find_if(inSegments.begin(), inSegments.end(), [sUId](InSegment* obj) {return obj->getuId() == sUId;}); // given a node uId, find it

    if (sId != inSegments.end()) { // gives us the segment index
        delete *sId;
        sIdx = std::distance(inSegments.begin(), sId);
        deleted[sIdx] = true;
        inSegments.erase(sId);

    }else{
        fprintf(stderr, "Error: cannot detect segment to be deleted (sHeader: %s). Terminating.\n", sHeader.c_str());
        exit(1);
    }
    
    return true;
    
}

void InSequences::removePath(unsigned int pUId, bool all, bool silent) {
    
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
    
    if (all) {
        
        std::vector<PathComponent> pathComponents = pathIt->getComponents();
        
        for (auto &component : pathComponents) {
            
            unsigned int cUId = component.id;
            
            if (component.componentType == SEGMENT) {

                auto sId = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                if (sId != inSegments.end())
                inSegments.erase(sId);

            }else if (component.componentType == GAP) {
                
                auto gId = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                if (gId != inGaps.end())
                inGaps.erase(gId);
                
            }
            
        }
        
    }
    
    if (pathIt != inPaths.end()) {
        
        inPaths.erase(pathIt);
        
    }else if (!silent) {
        
        fprintf(stderr, "Warning: the path you are attempting to remove does not exist (pUId: %i). Skipping.\n", pUId);
        
    }

}

void InSequences::removeGap(unsigned int gUId, bool silent) {
    
    auto gapIt = find_if(inGaps.begin(), inGaps.end(), [gUId](InGap& obj) {return obj.getuId() == gUId;}); // given a path pUId, find it
    
    if (gapIt != inGaps.end()) {
        
        inGaps.erase(gapIt);
        
    }else if (!silent){
        
        fprintf(stderr, "Warning: the gap you are attempting to remove does not exist (pUId: %i). Skipping.\n", gUId);
        
    }

}

void InSequences::resizeGap(std::string gHeader, unsigned int size) {
    
    auto gapIt = find_if(inGaps.begin(), inGaps.end(), [gHeader](InGap& obj) {return obj.getgHeader() == gHeader;}); // given a path pUId, find it
    
    if (gapIt != inGaps.end()) {
        
        gapIt->setDist(size);
        
    }else{
        
        fprintf(stderr, "Warning: the gap you are attempting to resize does not exist (pUId: %s). Skipping.\n", gHeader.c_str());
        
    }

}

void InSequences::removePathsFromSegment(unsigned int uId) {
    
    for (InPath& inPath : inPaths) {
        
        std::vector<PathComponent> pathComponents = inPath.getComponents();
        
        auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [uId](PathComponent& obj) {return obj.id == uId;}); // given a node uId, find if present in the given path
        
        if (pathIt != pathComponents.end()) {
        
            removePath(inPath.getpUId());

        }
        
    }
    
}

void InSequences::removePathComponents(unsigned int uId) {
    
    for (InPath& inPath : inPaths) {
        
        std::vector<PathComponent> pathComponents = inPath.getComponents();
        
        auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [uId](PathComponent& obj) {return obj.id == uId;}); // given a node uId, find if present in the given path
        
        if (pathIt != pathComponents.end()) {
        
            removePath(inPath.getpUId());

        }
        
    }
    
}

void InSequences::removeSegmentInPath(unsigned int suId, InGap gap) {
    
    std::vector<PathComponent> newComponents;
    
    for (InPath& inPath : inPaths) {
        
        std::vector<PathComponent> pathComponents = inPath.getComponents();
        
        auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [suId](PathComponent& obj) {return obj.id == suId;}); // given a node uId, find if present in the given path
        
        if (pathIt != pathComponents.end()) {
            
            newComponents.insert(newComponents.end(), std::begin(pathComponents), pathIt-1); // subset path excluding segment to be removed
            
            newComponents.push_back({GAP, gap.getuId(), '0', 0, 0}); // introduce gap
            
            newComponents.insert(newComponents.end(), pathIt+2, std::end(pathComponents)); // add remaining components
            
            inPath.setComponents(newComponents);
            
            break;
                
        }
        
    }
    
}

void InSequences::joinPaths(std::string pHeader, unsigned int pUId1, unsigned int pUId2, std::string gHeader, unsigned int gUId, char pId1Or, char pId2Or, unsigned int dist, unsigned int start1, unsigned int end1, unsigned int start2, unsigned int end2) {
    
    InPath path;
    
    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find(pHeader); // get the headers to uIds table to look for the header
    
    if (got == headersToIds.end()) { // this is the first time we see this path
        
        lg.verbose("Path not found in keys. Creating new path (" + pHeader + ", pUId: " + std::to_string(uId.get()) + ")");
        
        insertHash(pHeader, uId.get());
        
        path.newPath(uId.get(), pHeader);
        
        uId.next();
        
    }else{
        
        path.setHeader(pHeader);
        pUId1 = got->second;
        path.setpUId(pUId1);
        
        lg.verbose("Path already exists in keys (" + pHeader + ", pUId: " + std::to_string(pUId1) + "). Joining.");
        
    }
    
    PathComponent component1, component2;
        
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId1](InPath& obj) {return obj.getpUId() == pUId1;}); // given a path pUId, find it
    
    if (pathIt != inPaths.end()) {
        
        lg.verbose("Path found in path set (pIUd: " + std::to_string(pUId1) + "). Adding components to new path.");
        
        std::vector<PathComponent> pathComponents = pathIt->getComponents();
        
        trimPathByRef(pathComponents, start1, end1);
        
        if (pId1Or == '-') {revComPathComponents(pathComponents);}
        
        path.append({std::begin(pathComponents), std::end(pathComponents)});
        
        component1 = *std::prev(std::end(pathComponents));
        
        if (start1 == 0) {
        
            removePath(pUId1); // remove path1
            
        }
        
    }else{

        fprintf(stderr, "Error: could not locate in path set (pIUd: %u)\n", pUId1);
        exit(EXIT_FAILURE);
        
    }
    
    if (gUId == 0) {
        
        insertHash(gHeader, uId.get());
        
        gUId = uId.get();
        
        uId.next();
    
    }
    
    lg.verbose("Adding gap to new path (" + gHeader + ")");
    
    path.add(GAP, gUId, '0');
        
    pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId2](InPath& obj) {return obj.getpUId() == pUId2;}); // given a path pUId, find it
    
    if (pathIt != inPaths.end()) {
        
        lg.verbose("Path found in path set (pIUd: " + std::to_string(pUId2) + "). Adding components to new path.");
        
        std::vector<PathComponent> pathComponents = pathIt->getComponents();
        
        trimPathByRef(pathComponents, start2, end2);
        
        if (pId2Or == '-') {revComPathComponents(pathComponents);}
    
        path.append({std::begin(pathComponents), std::end(pathComponents)});
        
        component2 = *std::begin(pathComponents);
        
        if (start2 == 0) {
        
            removePath(pUId2); // remove path2
            
        }
        
    }else{
        
        fprintf(stderr, "Error: could not locate in path set (pIUd: %u)\n", pUId2);
        exit(EXIT_FAILURE);
        
    }
    
    InGap gap;
    
    Tag tag;
    memcpy(tag.label, "SC", sizeof tag.label);
    tag.type = 'i';
    tag.content = "1";
    
    gap.newGap(gUId, component1.id, component2.id, component1.orientation, component2.orientation, dist, gHeader, {tag}); // define the new gap
    
    addGap(gap); // introduce the new gap
    
    addPath(path);
    
}

InPath InSequences::joinPathsByComponent(std::string seqHeader, unsigned int uId1, unsigned int uId2, unsigned int uId3) {
    
    InPath newPath;
    newPath.setHeader(seqHeader);
    
    for (InPath inPath : inPaths) { // add first path
        
        std::vector<PathComponent> pathComponents = inPath.getComponents();
        
        auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [uId1](PathComponent& obj) {return obj.id == uId1;}); // given a node uId, find if present in the given path
        
        if (pathIt != pathComponents.end()) {
        
            newPath.append({std::begin(pathComponents), std::end(pathComponents)});
        
            break;
            
        }
        
    }
    
    newPath.add(GAP, uId2, '0');
    
    for (InPath inPath : inPaths) { // add second path
        
        std::vector<PathComponent> pathComponents = inPath.getComponents();
        
        auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [uId3](PathComponent& obj) {return obj.id == uId3;}); // given a node uId, find if present in the given path
        
        if (pathIt != pathComponents.end()) {
        
            newPath.append({std::begin(pathComponents), std::end(pathComponents)});
        
            break;
            
        }
        
    }
    
    return newPath;
    
}

void InSequences::splitPath(unsigned int guId, std::string pHeader1, std::string pHeader2) {
    
    InPath newPath1;
    newPath1.setHeader(pHeader1);
    std::vector<PathComponent> newComponents1;
    
    insertHash(pHeader1, uId.get());
    
    path.newPath(uId.get(), pHeader1);
    
    uId.next();
    
    InPath newPath2;
    newPath2.setHeader(pHeader2);
    std::vector<PathComponent> newComponents2;
    
    insertHash(pHeader2, uId.get());
    
    path.newPath(uId.get(), pHeader2);
    
    uId.next();
    
    for (InPath& inPath : inPaths) { // search through all paths
        
        std::vector<PathComponent> pathComponents = inPath.getComponents();
        
        auto pathIt = find_if(pathComponents.begin(), pathComponents.end(), [guId](PathComponent& obj) {return obj.id == guId;}); // given a node uId, find if present in the given path
        
        if (pathIt != pathComponents.end()) {
        
            newComponents1.insert(std::end(newComponents1), std::begin(pathComponents), pathIt);
            
            newPath1.setComponents(newComponents1);
            
            if (newComponents1.size() > 0) {
            
                addPath(newPath1);
                
            }
            
            newComponents2.insert(std::end(newComponents2), pathIt+1, std::end(pathComponents));
            
            newPath2.setComponents(newComponents2);
            
            if (newComponents2.size() > 0) {
            
                addPath(newPath2);
                
            }
            
            removePath(inPath.getpUId());
        
            break;
            
        }
        
    }
    
}

void InSequences::clearPaths() {
    
    inPaths.clear();
    
}

void InSequences::clearGaps() {
    
    inGaps.clear();
    
}

void InSequences::renamePath(unsigned int pUId, std::string pHeader, unsigned int* newpUId) {
    
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
    
    pathIt->setHeader(pHeader);
    
    insertHash(pHeader, pUId);
    
    if (newpUId != NULL) {
        
        pathIt->setpUId(*newpUId);
        
    }
    
}

void InSequences::revComPath(unsigned int pUId) {
    
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
    
    pathIt->revCom();
    
    lg.verbose("Path reverse-complemented (" + pathIt->getHeader() + ").");
    
}

void InSequences::trimPathByUId(unsigned int pUId, uint64_t start, uint64_t end) {
    
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
    
    std::vector<PathComponent>* pathComponents = pathIt->getComponentsByRef();
    
    trimPath(pathComponents, start, end);
    
}

void InSequences::trimPathByRef(std::vector<PathComponent>& pathComponents, uint64_t start, uint64_t end) {
    
    trimPath(&pathComponents, start, end);
    
}

void InSequences::trimPath(std::vector<PathComponent>* pathComponents, uint64_t start, uint64_t end) {
    
    if(start == 0 || end == 0) {

        lg.verbose("Nothing to trim. Skipping.");
        return;

    }
    
    lg.verbose("Trimming path (start: " + std::to_string(start) + ", end: " + std::to_string(end) + ")");
    
    std::string trimmed;
    
    uint64_t traversedSize = 0, actualSize = 0, compSize = 0, newCompSize = 0, compOriginalSize = 0;
    
    for (std::vector<PathComponent>::iterator component = pathComponents->begin(); component != pathComponents->end(); component++) {
        
        lg.verbose("New path size before iteration: " + std::to_string(actualSize));
        
        lg.verbose("Checking original coordinates of component (uId: " + std::to_string(component->id) + ", start: " + std::to_string(component->start) + ", end: " + std::to_string(component->end) + ")");
        
        compOriginalSize = getComponentSize(*component, true);
        
        compSize = getComponentSize(*component, false);
        
        trimmed = compSize != compOriginalSize ? " " : " not ";
            
        lg.verbose("Component was" + trimmed + "already trimmed (size: " + std::to_string(compSize) + ")");
            
        if (traversedSize + compSize < start) {
            
            lg.verbose("Start coordinate exceeds component, removing it");
            
            component--;
            
            pathComponents->erase(std::next(component));
            
            traversedSize += compSize;
            
            continue;
            
        }
		
		// this is where we could test if we are cutting in a segment, and if so, we could decide to avoid cutting inside the segment by realigning the cut
		// one option is to have a majority rule where the largest side takes all the sequence. Optionally there could also be a length cutoff for realignment
		// one challenge is that in some instances there is no gap, so a list of allowed cuts should be provided
        
        if (traversedSize + compSize >= start && traversedSize < start - 1 && traversedSize + compSize > end) {
           
           lg.verbose("Trimming both ends");
           
           trimComponent(*component, start - traversedSize, end - traversedSize);
           
        } else if (traversedSize + compSize >= start && traversedSize < start && traversedSize + compSize <= end) {
            
            lg.verbose("Trimming left end");
            
            trimComponent(*component, start - traversedSize, 0);
            
        } else if (traversedSize + compSize > end) {
            
            lg.verbose("Trimming right end");
            
            trimComponent(*component, 0, end - traversedSize);
            
        }
        
        if (traversedSize + compSize >= end) {
        
            pathComponents->erase(component + 1, pathComponents->end());
            
            lg.verbose("Erased extra components");
            
        }
        
        if (component->start > 0 && component->end == 0) { // account for editing of the start coordinate but not the end coordinate
            
            component->end = compOriginalSize;
            
        }else if (component->end > 0 && component->start == 0) { // account for editing of the end coordinate but not the start coordinate
            
            component->start = 1;
            
        }else if(component->start == 1 && component->end == compOriginalSize) { // if the result of trimming restores the original size of the component, no need to adjust the internal coordinates
            
            component->start = 0;
            component->end = 0;
            
        }
        
        newCompSize = getComponentSize(*component, false);
        lg.verbose("Component size after trimming: " + std::to_string(newCompSize));
        
        actualSize += newCompSize;
        lg.verbose("Path size after iteration: " + std::to_string(actualSize));
        
        traversedSize += compSize;
        lg.verbose("Traversed path: " + std::to_string(traversedSize));
        
    }
        
    lg.verbose("Final path size: " + std::to_string(actualSize));
    
    if (actualSize != end-start+1) {fprintf(stderr, "Error: Path size after trimming (%llu) differs from expected size after trimming (%llu). Terminating.\n", actualSize, end-start+1); exit(1);}

}

void InSequences::trimComponent(PathComponent& component, int start, int end) {
    
    int startCom = component.start, endCom = component.end;

    if (component.orientation == '+' || component.orientation == '0') { // we only change the end coordinate if the component wasn't already flipped, otherwise we edit the start
    
        if (start != 0 && end != 0) {
        
            component.end = (startCom == 0 ? 0 : startCom - 1) + end;
            
            component.start = (startCom == 0 ? 0 : startCom - 1) + start;
            
            lg.verbose("Plus orientation (+). Both start and end coordinates of the component need to be edited as result of subsetting (new start: " + std::to_string(component.start) + ", new end: " + std::to_string(component.end) + ")");
            
        }else if (start != 0) {
                    
            component.start = (startCom == 0 ? 0 : startCom - 1) + start;
                
            lg.verbose("Plus orientation (+). Start coordinate of the component needs to be edited as result of subsetting (new start: " + std::to_string(component.start) + ")");
 
        }else if (end != 0) {
         
            component.end = (startCom == 0 ? 0 : startCom - 1) + end;
            
            lg.verbose("Plus orientation (+). End coordinate of the component needs to be edited as result of subsetting (new end: " + std::to_string(component.end) + ")");
            
        }
        
    }else{
        
        if (start != 0 && end != 0) {
            
            component.start = (startCom == 0 ? 0 : startCom - 1) + getComponentSize(component, false) - end + 1;
            
            component.end = (endCom == 0 ? getComponentSize(component, false) : endCom) - start + 1;

            lg.verbose("Minus orientation (-). Both start and end coordinates of the component need to be edited as result of subsetting (new start: " + std::to_string(component.start) + ", new end: " + std::to_string(component.end) + ")");
            
        }else if (start != 0) {

            component.end = (endCom == 0 ? getComponentSize(component, false) : endCom) - start + 1;

            lg.verbose("Minus orientation (-). End coordinate of the component needs to be edited as result of subsetting (new end: " + std::to_string(component.end) + ")");

        }else if (end != 0) {

            component.start = (startCom == 0 ? 0 : startCom - 1) + getComponentSize(component, false) - end + 1;

            lg.verbose("Minus orientation (-). Start coordinate of the component needs to be edited as result of subsetting (new start: " + std::to_string(component.start) + ")");

        }

    }
    
}

uint64_t InSequences::getComponentSize(PathComponent& component, bool original) {
    
    unsigned int cUId = component.id;
    
    if (component.componentType == SEGMENT) {
    
        auto inSegment = find_if(inSegments.begin(), inSegments.end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
        
        return original ? (*inSegment)->getSegmentLen() : (*inSegment)->getSegmentLen(component.start, component.end);
        
    }else{
        
        auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
        
        return original ? inGap->getDist() : inGap->getDist(component.start, component.end);
        
    }
    
}

unsigned int InSequences::pathLen(unsigned int pUId) {
    
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
    
    return pathIt->getLen();
    
}

InSegment& InSequences::findSegmentBySUId(uint32_t sUId) {
    
    auto inSegment = find_if(inSegments.begin(), inSegments.end(), [sUId](InSegment* obj) {return obj->getuId() == sUId;}); // given a node Uid, find it
    
    if (inSegment != inSegments.end()) {
        return **inSegment;
    }else{
        fprintf(stderr, "Error: segment sUId not found (sUId: %u). Terminating.\n", sUId); exit(1);
    }
    
}

void InSequences::walkPath(InPath* path) {
    
    unsigned int cUId = 0, gapLen = 0;
    
    std::vector<PathComponent> pathComponents = path->getComponents();
    
    path->reinitializeCounts();
    
    for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
        cUId = component->id;
    
        if (component->componentType == SEGMENT) {
            
            InSegment& inSegment = findSegmentBySUId(cUId);
            
            contigLens.push_back(inSegment.getSegmentLen(component->start, component->end));
            
            path->increaseContigN();
            
            path->increaseLen(inSegment.getSegmentLen(component->start, component->end));
            
            path->increaseSegmentLen(inSegment.getSegmentLen(component->start, component->end));
            
            path->increaseLowerCount(inSegment.getLowerCount(component->end - component->start));
            
            if (component->start == 0 || component->end == 0) {
            
                path->increaseA(component->orientation == '+' ? inSegment.getA() : inSegment.getT());
                
                path->increaseC(component->orientation == '+' ? inSegment.getC() : inSegment.getG());
                
                path->increaseG(component->orientation == '+' ? inSegment.getG() : inSegment.getC());
                
                path->increaseT(component->orientation == '+' ? inSegment.getT() : inSegment.getA());
                
            }else{
                
                std::string sequence = inSegment.getInSequence(component->start, component->end);
                
                if (component->orientation == '+') {
                
                    for (char base : sequence) {
                        
                        switch (base) {
                            case 'A':
                            case 'a':{
                                
                                path->increaseA(1);
                                break;
                                
                            }
                            case 'C':
                            case 'c':{
                                
                                path->increaseC(1);
                                break;
                                
                            }
                            case 'G':
                            case 'g': {
                                
                                path->increaseG(1);
                                break;
                                
                            }
                            case 'T':
                            case 't': {
                                
                                path->increaseT(1);
                                break;
                                
                            }
                    
                        }
                        
                    }

                }else{
                    
                    for (char base : sequence) {
                    
                        switch (base) {
                            case 'A':
                            case 'a':{
                                
                                path->increaseT(1);
                                break;
                                
                            }
                            case 'C':
                            case 'c':{
                                
                                path->increaseG(1);
                                break;
                                
                            }
                            case 'G':
                            case 'g': {
                                
                                path->increaseC(1);
                                break;
                                
                            }
                            case 'T':
                            case 't': {
                                
                                path->increaseA(1);
                                break;
                                
                            }
                    
                        }
                        
                    }
                    
                }
                
            }
            
        }else if (component->componentType == GAP){
            
            auto inGap = find_if(inGaps.begin(), inGaps.end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
            
            gapLen += inGap->getDist(component->start, component->end);
            
            if (component + 1 == pathComponents.end() || !((component + 1)->componentType == GAP)) {
            
                gapLens.push_back(gapLen);
                
                gapLen = 0;
                
            }
            
            path->increaseLen(inGap->getDist(component->start, component->end));
            
        }else{} // need to handle edges, cigars etc
        
    }
    
}

void InSequences::discoverPaths() {
    
    buildGraph(inGaps);
    
    for (InSegment* inSegment : inSegments) {
        
        if (!visited[inSegment->getuId()]) {
            
            InPath path;
            
            path.newPath(uId.get(), inSegment->getSeqHeader() + "_path");
            
            insertHash(inSegment->getSeqHeader() + "_path", uId.get());
            
            uId.next();
            
            dfsPath(inSegment->getuId(), path);
            
            addPath(path);
            
        }
        
    }
    
}

unsigned int InSequences::detectOverlap(std::string* sequence1, std::string* sequence2, unsigned int minOvlLen) {
    
    unsigned int ovlLen = 0;
    
    if (sequence1->size() > minOvlLen && sequence2->size() > minOvlLen) {
        
        std::size_t found = 0, pos = 0;
        
        while (found < std::string::npos) {
            
            found = sequence2->find(sequence1->substr(sequence1->size()-minOvlLen, minOvlLen), pos);
            
            if (found == std::string::npos)
                break;
            
            if (sequence1->substr(sequence1->size()-minOvlLen-found, minOvlLen+found) == sequence2->substr(0, found+minOvlLen))
                ovlLen = minOvlLen+found;
                
            pos = minOvlLen+found;
            
        }
        
    }
    
    if (ovlLen == sequence1->size() || ovlLen == sequence2->size())
        return 0;
    
    return ovlLen;
    
}

void InSequences::discoverTerminalOverlaps(int terminalOvlLen) {
    
    unsigned int max = inSegments.size();
    InSegment* inSegment1, *inSegment2;
        
    for (unsigned int i = 0; i < max; ++i) {
        
        for (unsigned int j = i; j < max; ++j) {
            
            inSegment1 = inSegments[i];
            inSegment2 = inSegments[j];
            
            std::string* sequence1 = inSegment1->inSequence;
            std::string* sequence2 = inSegment2->inSequence;
            
            unsigned int maxOverhang = 20000;
            
            if(sequence1->size() < maxOverhang)
                maxOverhang = sequence1->size();
            
            if(sequence2->size() < maxOverhang && sequence2->size() < sequence1->size())
                maxOverhang = sequence2->size();
            
            std::string subSeq1 = sequence1->substr(sequence1->size()-maxOverhang, maxOverhang);
            std::string subSeq2 = sequence2->substr(0, maxOverhang);
            
            std::string subSeq3 = revCom(sequence2->substr(sequence2->size()-maxOverhang, maxOverhang));
            std::string subSeq4 = revCom(sequence1->substr(0, maxOverhang));
            
            unsigned int ovlLen = 0;
            
            ovlLen = detectOverlap(&subSeq1, &subSeq2, terminalOvlLen);
            
            if (ovlLen>0) {
                
                lg.verbose("Found perfect overlap of length " + std::to_string(ovlLen) + " between: " + inSegment1->getSeqHeader() + "+ and " + inSegment2->getSeqHeader() + "+");
                
                edge.newEdge(this->uId.next(), inSegment1->getuId(), inSegment2->getuId(), '+', '+', std::to_string(ovlLen) + "M");
                
                this->appendEdge(edge);
                
            }
            
            ovlLen = detectOverlap(&subSeq1, &subSeq3, terminalOvlLen);

            if (ovlLen>0) {
             
                lg.verbose("Found perfect overlap of length " + std::to_string(ovlLen) + " between: " + inSegment1->getSeqHeader() + "+ and " + inSegment2->getSeqHeader() + "-");
                
                edge.newEdge(this->uId.next(), inSegment1->getuId(), inSegment2->getuId(), '+', '-', std::to_string(ovlLen) + "M");
                
                this->appendEdge(edge);
                
            }

            ovlLen = detectOverlap(&subSeq4, &subSeq2, terminalOvlLen);

            if (ovlLen>0) {
                lg.verbose("Found perfect overlap of length " + std::to_string(ovlLen) + " between: " + inSegment1->getSeqHeader() + "- and " + inSegment2->getSeqHeader() + "+");
             
                edge.newEdge(this->uId.next(), inSegment1->getuId(), inSegment2->getuId(), '-', '+', std::to_string(ovlLen) + "M");
                
                this->appendEdge(edge);
                
            }
            
        }
        
    }
    
}

void InSequences::dfsPath(unsigned int v, InPath& newPath) // Depth First Search to build a new path given a vertex
{

    visited[v] = true; // mark the current node as visited

    if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(adjListFW.at(v).at(0).segmentId == adjListBW.at(v).at(0).segmentId) && !backward) { // if the vertex has exactly one forward and one backward connection and they do not connect to the same vertex (internal node)

        lg.verbose("node: " + idsToHeaders[v] + " --> case a: internal node, forward direction");
        
        newPath.add(SEGMENT, v, '+');

        backward = false;

    }else if (adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 1){ // this is the final vertex without gaps

        lg.verbose("node: " + idsToHeaders[v] + " --> case b: end node, forward direction, no final gap");
        
        newPath.add(SEGMENT, v, '-');

        backward = true; // reached the end

    }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 2){ // this is the final vertex with terminal gap

        lg.verbose("node: " + idsToHeaders[v] + " --> case c: end node, forward direction, final gap");
        
        newPath.add(SEGMENT, v, '-');

        if (adjListBW.at(v).at(0).segmentId != v) { // make sure you are not using the terminal edge to ascertain direction in case it was edited by sak

        }else{

        }

        backward = true; // reached the end

    }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && !(adjListFW.at(v).at(0).segmentId == adjListBW.at(v).at(0).segmentId) && backward){ // this is an intermediate vertex, only walking back

        lg.verbose("node: " + idsToHeaders[v] + " --> case d: intermediate node, backward direction");
        
        newPath.add(SEGMENT, v, '-');

        backward = true;

    }else if(adjListFW.at(v).size() == 0 && adjListBW.at(v).size() == 0){ // disconnected component
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case e: disconnected component");
        
        newPath.add(SEGMENT, v, '+');

    }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 0){ // this is the first vertex without gaps
        
        lg.verbose("node: " + idsToHeaders[v] + " --> case f: start node, no gaps");
        
        newPath.add(SEGMENT, v, '+');

        visited.clear();

        visited[v] = true; // we have just visited the start node

        backward = false;

    }else if (adjListFW.at(v).size() == 2 && adjListBW.at(v).size() == 1){ // this is the first vertex with a terminal gap

        lg.verbose("node: " + idsToHeaders[v] + " --> case g: start node, start gap");

        newPath.add(SEGMENT, v, '+');

        visited.clear();

        visited[v] = true; // we have just visited the start node

        backward = false;

    }else if (adjListFW.at(v).size() == 1 && adjListBW.at(v).size() == 1 && adjListFW.at(v).at(0).segmentId == adjListBW.at(v).at(0).segmentId) { // if the vertex has exactly one forward and one backward connection and they connect to the same vertex (disconnected component with gap)

        lg.verbose("node: " + idsToHeaders[v] + " --> case h: disconnected component with gap");
        
        newPath.add(SEGMENT, v, '+');

        backward = false;

    }

    for (auto i: adjListFW[v]) { // recur for all forward vertices adjacent to this vertex

        if (!visited[i.segmentId] && !deleted[i.segmentId]) {

            dfsPath(i.segmentId, newPath); // recurse

        }
    }

    for (auto i: adjListBW[v]) { // recur for all backward vertices adjacent to this vertex

        if (!visited[i.segmentId] && !deleted[i.segmentId]) {

            dfsPath(i.segmentId, newPath); // recurse

        }
    }

}

void InSequences::findBubbles() {
    
    unsigned int sUId = 0, sUId1 = 0, sUId2 = 0;
    char sId1Or, sId2Or;
    
    visited.clear();
    
    buildEdgeGraph();
    
    for (InSegment* segment : inSegments) {
        
        sUId = segment->getuId();
        
        lg.verbose("Evaluating for node: " + idsToHeaders[sUId] + " (uId: " + std::to_string(sUId) + ")");
        
        if (!visited[sUId] && !deleted[sUId]) { // check if the node was already visited
            
            lg.verbose("The node was not yet visited or deleted");
            
            if (adjEdgeList.at(sUId).size() == 2 && adjEdgeList.at(sUId).at(0).orientation0 != adjEdgeList.at(sUId).at(1).orientation0) { // if it has exactly two edges with different orientation it could be a het region of the bubble
                
                lg.verbose("Exactly two edges with different orientation found. Could be a het region of the bubble");
                
                // then check the the adjacient nodes
                sUId1 = adjEdgeList.at(sUId).at(0).id;
                sId1Or = adjEdgeList.at(sUId).at(0).orientation1;
                
                sUId2 = adjEdgeList.at(sUId).at(1).id;
                sId2Or = adjEdgeList.at(sUId).at(1).orientation1;
                
                if (adjEdgeList.at(sUId1).size() >= 2 && adjEdgeList.at(sUId2).size() >= 2) { // both nodes need at least two edges to be a bubble
                    
                    lg.verbose("Both neighbour nodes have at least two edges");
                    
                    for (auto edge1 : adjEdgeList.at(sUId1)) {
                        
                        if (edge1.orientation1 == sId1Or && edge1.id != sUId) { // we are checking edges on the side of the potential bubble for node1, avoiding the original node
                            
                            lg.verbose("Evaluating node: " + idsToHeaders[edge1.id] + " (uId: " + std::to_string(edge1.id) + ")");
                            
                            if (edge1.id == sUId2) { // this is a potential insertion
                                
                                bubbles.push_back({sUId, sUId1, sUId2, 0});
                                
                                lg.verbose("Candidate insertion found");
                                
                            }else{
                                
                                for (auto edge2 : adjEdgeList.at(sUId2)) {
                                
                                    if (edge2.orientation1 == sId2Or && edge1.id == edge2.id) { // we are checking edges on the side of the potential bubble for node2 and that it connects to the same node as node1
                                        
                                        lg.verbose("Evaluating node: " + idsToHeaders[edge2.id] + " (uId: " + std::to_string(edge2.id) + ")");
                                        
                                        bubbles.push_back({sUId, sUId1, sUId2, edge2.id});
                                        
                                        lg.verbose("Candidate bubble found");
                                        
                                    }
                                    
                                }
                                
                            }
                               
                        }
                        
                    }
                    
                }
                
            }

        }
        
        visited[sUId] = true;
        
    }
    
}

std::vector<Bubble>* InSequences::getBubbles() {
    return &bubbles;
}

std::vector<uint32_t> InSequences::getCircularSegments() {
    
    std::vector<uint32_t> segments;
    std::vector<uint32_t>::iterator it;
    
    for (InEdge& inEdge : inEdges) {
        if (inEdge.sId1 == inEdge.sId2)
            segments.push_back(inEdge.sId1);
    }
    sort(segments.begin(), segments.end());
    it = std::unique(segments.begin(), segments.end());
    segments.resize(std::distance(segments.begin(),it));
    return segments;
}

std::vector<uint32_t> InSequences::getCircularPaths() {
    
    std::vector<uint32_t> pathIds;
    std::vector<uint32_t>::iterator it;
    std::vector<uint32_t> circularSegmentsIds = getCircularSegments();
    
    for (InPath& inPath : inPaths) {
        
        std::vector<PathComponent> pathComponents = inPath.getComponents();
        
        if (pathComponents.size() > 0 && (pathComponents.front().id == pathComponents.back().id)) {
            
            uint32_t componentId = pathComponents.begin()->id;

            if (std::find(circularSegmentsIds.begin(), circularSegmentsIds.end(), componentId) != circularSegmentsIds.end())
                pathIds.push_back(inPath.pUId);
            

            
        }
        
    }
    
    sort(pathIds.begin(), pathIds.end());
    it = std::unique(pathIds.begin(), pathIds.end());
    pathIds.resize(std::distance(pathIds.begin(),it));
    
    return pathIds;
    
}

void InSequences::maskPath(std::string pHeader, unsigned int start, unsigned int end, unsigned int dist) {
    
    unsigned int pUId = 0;
    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find(pHeader); // get the headers to uIds table to look for the header
    
    if (got == headersToIds.end()) { // this is the first time we see this path
        fprintf(stderr, "Error: path name not found (%s). Terminating.\n", pHeader.c_str()); exit(1);
    }else{
        pUId = got->second;
    }
    
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;}); // given a path pUId, find it
    std::vector<PathComponent>* pathComponents = pathIt->getComponentsByRef();
    
    if(start == 0 || end == 0) {
        lg.verbose("Nothing to mask. Skipping.");
        return;
    }
    
    lg.verbose("Masking path (start: " + std::to_string(start) + ", end: " + std::to_string(end) + ")");
    
    uint64_t traversedSize = 0, actualSize = 0, compSize = 0, newCompSize = 0;
    
    for (std::vector<PathComponent>::iterator component = pathComponents->begin(); component != pathComponents->end(); component++) {
        
        walkPath(&(*pathIt));
        
        lg.verbose("Path size before iteration: " + std::to_string(pathIt->getLen()));
        
        lg.verbose("Checking original coordinates of component (uId: " + std::to_string(component->id) + ", start: " + std::to_string(component->start) + ", end: " + std::to_string(component->end) + ")");
        
        compSize = getComponentSize(*component, false);
        
        if (traversedSize + compSize > start) {
            
            if (traversedSize + compSize < end) {fprintf(stderr, "Error: end coordinate exceeds length of component (%u). Terminating.\n", end); exit(1);}
            
            lg.verbose("Component to mask identified. Splitting");
            
            component->end = traversedSize + component->start + start;
            
            component->start = component->start + 1;
            
            lg.verbose("Adding gap of size: " + std::to_string(dist));
            
            InGap gap;
            
            gap.newGap(uId.get(), component->id, component->id, '+', '+', dist, getInSegment(component->id)->getSeqHeader() + ".innerGap");
            
            addGap(gap); // introduce the new gap
            
            unsigned int gUId = gap.getuId(), sUId = component->id, newStart = traversedSize + component->start + end, newEnd = compSize;
            char orientation = component->orientation;
            
            insertHash(gap.getgHeader(), gUId);
            
            pathComponents->insert(std::next(component), {GAP, gUId, '0', 0, 0});
            
            lg.verbose("Adding residual component of size: " + std::to_string(newEnd - newStart));
            
            auto componentIt = find_if(pathComponents->begin(), pathComponents->end(), [gUId](PathComponent& obj) {return obj.id == gUId;});
            
            pathComponents->insert(std::next(componentIt), {SEGMENT, sUId, orientation, newStart, newEnd});
        
            newCompSize = getComponentSize(*component, false);
            lg.verbose("Component size before gap: " + std::to_string(newCompSize));
            lg.verbose("Gap component size: " + std::to_string(dist));
            lg.verbose("Component size after gap: " + std::to_string(newEnd - newStart + 1));
            
            walkPath(&(*pathIt));
            lg.verbose("Path size after iteration: " + std::to_string(pathIt->getLen()));

            break;
            
        }
        
    }
        
    lg.verbose("Final path size: " + std::to_string(actualSize));
    
    //if (actualSize != end-start+1) {fprintf(stderr, "Error: Path size after trimming (%u) differs from expected size after trimming (%u). Terminating.\n", actualSize, end-start+1); exit(1);}

}

void InSequences::updateEdgeSUId(uint32_t sUId, uint32_t new_sUId, char vertex) {
    
    for (InEdge &edge: inEdges) {
        
        if (vertex == 'L') {
            
            if (edge.getsId1() == sUId && edge.getsId1Or() == '-')
                edge.setsId1(new_sUId);
            
            if (edge.getsId2() == sUId && edge.getsId1Or() == '+')
                edge.setsId2(new_sUId);
            
        }else if (vertex == 'R') {
            
            if (edge.getsId1() == sUId && edge.getsId1Or() == '+')
                edge.setsId1(new_sUId);
            
            if (edge.getsId2() == sUId && edge.getsId1Or() == '-')
                edge.setsId2(new_sUId);
            
        }else{
            
            fprintf(stderr, "Error: unknown vertex (%c). Terminating.\n", vertex); exit(1);
            
        }
        
    }
    
}

std::pair<InSegment*,InSegment*> InSequences::cleaveSegment(uint32_t sUId, uint64_t start, std::string sHeader2, std::string sHeader3, std::string eHeader1) {
    
    InSegment *inSegment1 = getInSegment(sUId);
    InSegment *inSegment2 = new InSegment(*inSegment1);
    
    lg.verbose("Segment1 size before trimming: " + std::to_string(inSegment1->getSegmentLen()));
    inSegment1->trimSegment(start, inSegment1->getSegmentLen());
    inSegment1->setSeqHeader(sHeader2);
    inSegment1->setuId(uId.get());
    updateEdgeSUId(sUId, uId.get(), 'L');
    insertHash(sHeader2, uId.get());
    uId.next();
    lg.verbose("Segment1 size after trimming: " + std::to_string(inSegment1->getSegmentLen()));

    lg.verbose("Segment2 size before trimming: " + std::to_string(inSegment2->getSegmentLen()));
    inSegment2->trimSegment(0, start);
    inSegment2->setSeqHeader(sHeader3);
    inSegment2->setuId(uId.get());
    updateEdgeSUId(sUId, uId.get(), 'R');
    insertHash(sHeader3, uId.get());
    uId.next();
    inSegments.push_back(inSegment2);
    lg.verbose("Segment2 size after trimming: " + std::to_string(inSegment2->getSegmentLen()));
    
    if (eHeader1 != "") {
        edge.newEdge(this->uId.next(), inSegment1->getuId(), inSegment2->getuId(), '+', '+', "0M", eHeader1);
        this->appendEdge(edge);
    }
    return std::make_pair(inSegment1, inSegment2);
}

void InSequences::pushBackSegment(InSegment *inSegment) {
    InSegment* inSegmentCpy = new InSegment(*inSegment);
    inSegmentCpy->setuId(uId.next());
    inSegments.push_back(inSegmentCpy);
    insertHash(inSegmentCpy->seqHeader, inSegmentCpy->getuId());
}

InSequences* InSequences::subgraph(std::vector<std::string> nodeList) {
    
    phmap::flat_hash_set<std::string> nodes(nodeList.begin(), nodeList.end());
    
    InSequences *subgraph = new InSequences;
    
    for (InSegment *inSegment : inSegments) {
        if (nodes.find(inSegment->seqHeader) != nodes.end()) {
            subgraph->pushBackSegment(inSegment);
            subgraph->insertHash(inSegment->seqHeader, subgraph->uId.next());
        }
    }
    for (const InEdge &inEdge : inEdges) {
        
        std::string sHeader1 = idsToHeaders.find(inEdge.sId1)->second, sHeader2 = idsToHeaders.find(inEdge.sId2)->second;
        if (nodes.find(sHeader1) != nodes.end() && nodes.find(sHeader2) != nodes.end()) {
            InEdge inEdgeCpy(inEdge);
            inEdgeCpy.seteUId(subgraph->uId.next());
            inEdgeCpy.setsId1(subgraph->headersToIds.find(sHeader1)->second);
            inEdgeCpy.setsId2(subgraph->headersToIds.find(sHeader2)->second);
            subgraph->insertHash(inEdgeCpy.eHeader, inEdgeCpy.geteUId());
            subgraph->appendEdge(inEdgeCpy);
        }
    }
    return subgraph;
}

void InSequences::renamePath(std::string pHeader, std::string newHeader) {
    
    unsigned int pUId = 0;
    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find(pHeader); // get the headers to uIds table to look for the header
    
    if (got == headersToIds.end()) { // this is the first time we see this path
        fprintf(stderr, "Error: path name not found (%s). Terminating.\n", pHeader.c_str());
        exit(1);
    }else{
        pUId = got->second;
    }
    
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;});
    pathIt->setHeader(newHeader);
    eraseHash(pHeader, pUId);
    insertHash(newHeader, pUId);
}

void InSequences::updateComment(std::string pHeader, std::string comment) {
    
    unsigned int pUId = 0;
    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find(pHeader); // get the headers to uIds table to look for the header
    
    if (got == headersToIds.end()) { // this is the first time we see this path
        fprintf(stderr, "Error: path name not found (%s). Terminating.\n", pHeader.c_str());
        exit(1);
    }else{
        pUId = got->second;
    }
    auto pathIt = find_if(inPaths.begin(), inPaths.end(), [pUId](InPath& obj) {return obj.getpUId() == pUId;});
    pathIt->setComment(comment);
}
