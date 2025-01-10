#include <stdlib.h>
#include <string>
#include <vector>

#include <parallel-hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "log.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "stream-obj.h"
#include "gfa.h"
#include "input-filters.h"

Sequence* includeExcludeSeq(std::string seqHeader, std::string seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates bedExcludeList, std::string* inSequenceQuality) {
    
    uint64_t cBegin = 0, cEnd = 0, offset = 0, prevCEnd = 0;
    bool outSeq = false;
    
    lg.verbose("Processing sequence: " + seqHeader);
    
    if(!bedIncludeList.empty() || !bedExcludeList.empty())
        inSequence->erase(remove(inSequence->begin(), inSequence->end(), '\n'), inSequence->end());
    
    if   (bedIncludeList.empty() &&
          bedExcludeList.empty()) {
        
        outSeq = true;
        
    }else if(!bedIncludeList.empty() &&
              bedExcludeList.empty()) {
        
        offset = 0, prevCEnd = 0;
        outSeq = false;
        
        auto coordinates = bedIncludeList.getCoordinates();
        auto got = coordinates.find(seqHeader);
        if (got == coordinates.end())
            return NULL;
        
        for(std::pair<uint64_t,uint64_t> coordinate : got->second) {
            
            outSeq = true;
            cBegin = coordinate.first;
            cEnd = coordinate.second;
            
            if (!(cBegin == 0 && cEnd == 0)) {
                
                inSequence->erase(offset, cBegin-prevCEnd);
                
                if (inSequenceQuality != NULL)
                    inSequenceQuality->erase(offset, cBegin-prevCEnd);
                    
                offset += cEnd-cBegin;
                prevCEnd = cEnd;
                
            }
        }
            
        if (outSeq && inSequence->size()>0) {
            
            if (offset>0) {
            
                inSequence->erase(offset, inSequence->size()-offset);
                
                if (inSequenceQuality != NULL)
                    inSequenceQuality->erase(offset, inSequenceQuality->size()-offset);
            }
            outSeq = true;
        }
    }else if(bedIncludeList.empty() &&
            !bedExcludeList.empty()) {
        
        offset = 0;
        outSeq = true;
        
        auto coordinates = bedExcludeList.getCoordinates();
        auto got = coordinates.find(seqHeader);
        if (got == coordinates.end())
            return NULL;
        
        for(std::pair<uint64_t,uint64_t> coordinate : got->second) {

            cBegin = coordinate.first;
            cEnd = coordinate.second;
            
            if (!(cBegin == 0 && cEnd == 0)) {
                
                inSequence->erase(cBegin-offset, cEnd-cBegin);
                
                if (inSequenceQuality != NULL)
                    inSequenceQuality->erase(cBegin-offset, cEnd-cBegin);
                    
                offset += cEnd-cBegin;
                
            }else{
                outSeq = false;
            }
        }
                
    }else if
            (!bedIncludeList.empty() &&
             !bedExcludeList.empty()) {
                
                auto coordinates = bedIncludeList.getCoordinates();
                auto got1 = coordinates.find(seqHeader);
                coordinates = bedExcludeList.getCoordinates();
                auto got2 = coordinates.find(seqHeader);
                if (got1 == coordinates.end() && got2 == coordinates.end()) {
                    outSeq = true;
                }
    }
    if (outSeq && inSequence->size()>0) {
        return new Sequence{seqHeader, seqComment, *inSequence, *inSequenceQuality};
    }else {
        lg.verbose("Sequence entirely removed as a result of BED filter: " + seqHeader);
        return NULL;
    }
}

Sequence* includeExcludeSeg(InSequences* inSequences, std::string* seqHeader, std::string* seqComment, std::string* inSequence, BedCoordinates bedIncludeList, BedCoordinates* bedExcludeList, std::string* inSequenceQuality) {
    
    std::vector<std::string> bedIncludeListHeaders;
    std::vector<std::string> bedExcludeListHeaders;
    uint64_t cBegin = 0, cEnd = 0, offset = 0, prevCEnd = 0;

    bedIncludeListHeaders = bedIncludeList.getHeaders();
    
    if(bedExcludeList != NULL) {
    
        bedExcludeListHeaders = bedExcludeList->getHeaders();
    
    }
        
    bool outSeq = false;
    
    lg.verbose("Processing sequence: " + *seqHeader);
    
    if   (bedIncludeList.empty() &&
          (bedExcludeList == NULL ||
          bedExcludeList->empty())) {
        
        outSeq = true;
        
    }else if(!bedIncludeList.empty() &&
              bedExcludeList != NULL &&
              bedExcludeList->empty()) {
        
        if(inSequences->getInSegments()->size() == bedIncludeList.size()) { // check if we retrieved all we needed
            
            lg.verbose("Found all sequences, stop streaming input");
            outSeq = true;
            
        }
        
        offset = 0, prevCEnd = 0;
        outSeq = false;
        
        auto coordinates = bedIncludeList.getCoordinates();
        auto got = coordinates.find(*seqHeader);
        if (got == coordinates.end())
            return NULL;
        
        for(std::pair<uint64_t,uint64_t> coordinate : got->second) {
            
            outSeq = true;
            cBegin = coordinate.first;
            cEnd = coordinate.second;
            
            if (!(cBegin == 0 && cEnd == 0)) {
                
                inSequence->erase(offset, cBegin-prevCEnd);
                
                if (inSequenceQuality != NULL)
                    inSequenceQuality->erase(offset, cBegin-prevCEnd);
                    
                offset += cEnd-cBegin;
                prevCEnd = cEnd;
            }
        }
            
        if (outSeq && inSequence->size()>0) {
            
            if (offset>0) {
            
                inSequence->erase(offset, inSequence->size()-offset);
                
                if (inSequenceQuality != NULL)
                    inSequenceQuality->erase(offset, inSequenceQuality->size()-offset);
            }
            outSeq = true;
        }else{
            lg.verbose("Sequence entirely removed as a result of include: " + *seqHeader);
        }
    }else if(bedIncludeList.empty() &&
             bedExcludeList != NULL &&
            !bedExcludeList->empty()) {
            
        offset = 0;
        outSeq = true;
        
        auto coordinates = bedExcludeList->getCoordinates();
        auto got = coordinates.find(*seqHeader);
        if (got == coordinates.end())
            return NULL;
        
        for(std::pair<uint64_t,uint64_t> coordinate : got->second) {
            
            cBegin = coordinate.first;
            cEnd = coordinate.second;
            
            if (!(cBegin == 0 && cEnd == 0)) {
                
                inSequence->erase(cBegin-offset, cEnd-cBegin);
                
                if (inSequenceQuality != NULL)
                    inSequenceQuality->erase(cBegin-offset, cEnd-cBegin);
                    
                offset += cEnd-cBegin;
            }else{
                outSeq = false;
            }
        }
            
        if (outSeq && inSequence->size()>0)
            outSeq = true;
        else
            lg.verbose("Sequence entirely removed as a result of exclude: " + *seqHeader);
                
    }else if
            (!bedIncludeList.empty() &&
              bedExcludeList != NULL &&
             !bedExcludeList->empty()) {
                
                auto coordinates = bedIncludeList.getCoordinates();
                auto got1 = coordinates.find(*seqHeader);
                coordinates = bedExcludeList->getCoordinates();
                auto got2 = coordinates.find(*seqHeader);
                if (got1 == coordinates.end() && got2 == coordinates.end()) {
                    if(inSequences->getInSegments()->size() == bedIncludeList.size()) // check if we retrieved all we needed
                        lg.verbose("Found all sequences, stop streaming input");
                }
    }
    
    if (outSeq && inSequence->size()>0) {
        return new Sequence {*seqHeader, seqComment != NULL ? *seqComment : "", *inSequence, std::string()};
    }else {
        lg.verbose("Sequence entirely removed as a result of BED filter: " + *seqHeader);
        return NULL;
    }
}
