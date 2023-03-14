#ifndef KMER_H
#define KMER_H

#include <parallel_hashmap/phmap.h>
#include "parallel_hashmap/phmap_dump.h"

#include <fastx.h>

template<typename T>
inline void freeContainer(T& p_container)
{
    T empty;
    using std::swap;
    swap(p_container, empty);
}

template<typename VALUE>
struct Buf {
    uint64_t pos = 0, size = 10000;
    VALUE *seq = new VALUE[size];
};

template<class INPUT, typename VALUE, typename TYPE>
class Kmap {

protected:
    
    InSequences inSequences;
    
    uint32_t processedBuffers = 0;

    uint8_t k;
    
    uint64_t totKmers = 0, totKmersUnique = 0, totKmersDistinct = 0;
    
    const uint64_t mapCount = k < 28 ? pow(4,k/4) : pow(4,6);
    
    const uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount;
    
    uint64_t* pows = new uint64_t[k];

    std::vector<Buf<TYPE>*> buffers;
    
    phmap::flat_hash_map<uint64_t, VALUE>* map = new phmap::flat_hash_map<uint64_t, VALUE>[mapCount];
    
    std::vector<bool> mapsInUse = std::vector<bool>(mapCount, false);
    
    phmap::flat_hash_map<uint64_t, uint64_t> finalHistogram;
    
    const uint8_t ctoi[256] = {
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
    
    const uint8_t itoc[4] = {'A', 'C', 'G', 'T'};
    
public:
    
    std::vector<Log> logs;
    
    Kmap(uint8_t k) : k(k) {
        
        for(uint8_t p = 0; p<k; ++p)
            pows[p] = (uint64_t) pow(4,p);
        
    };
    
    ~Kmap(){
        
        delete[] map;
        delete[] pows;
        
    }
    
    void load(INPUT& userInput);
    
    void kunion(INPUT& userInput);
    
    bool mergeMaps(std::vector<std::string> prefixes, uint16_t m);
    
    bool unionSum(phmap::flat_hash_map<uint64_t, VALUE>& map1, phmap::flat_hash_map<uint64_t, VALUE>& map2);
    
    bool traverseInReads(Sequences* readBatch);
    
    bool traverseInReads(std::string* readBatch);
    
    void appendSequence(Sequence* sequence);
    
    void appendReads(Sequences* readBatch);
    
    inline uint64_t hash(uint8_t* string, bool* isFw = NULL);
    
    void hashSequences(Sequences* readBatch);
    
    void hashSequences(std::string* readBatch);
    
    void hashSegments();
    
    void consolidate();
    
    void finalize();
    
    void hist();
    
    bool countBuff(Buf<uint64_t>* buf, uint16_t m);
    
    bool countBuffs(uint16_t m);
    
    bool histogram(phmap::flat_hash_map<uint64_t, VALUE>& map);
    
    void resizeBuff(Buf<uint64_t>* buff);
    
    void printHist(std::unique_ptr<std::ostream>& ostream);
    
    void report(INPUT& userInput);
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);

};

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::appendSequence(Sequence* sequence) { // method to append a new sequence from a fasta
    
    threadPool.queueJob([=]{ return inSequences.traverseInSequence(sequence); });
    
    if(verbose_flag) {std::cerr<<"\n";};
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    lck.unlock();
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::load(INPUT& userInput){
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return loadMap(userInput.iSeqFileArg, m); });
    
    jobWait(threadPool);
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::kunion(INPUT& userInput){
        
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return mergeMaps(userInput.iReadFileArg, m); });
    
    jobWait(threadPool);
    
}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::mergeMaps(std::vector<std::string> prefixes, uint16_t m) {
    
    std::string prefix = prefixes[0]; // load the first map
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    map[m].phmap_load(ar_in);
    
    unsigned int numFiles = prefixes.size();

    for (unsigned int i = 1; i < numFiles; i++) {
        
        std::string prefix = prefixes[i]; // load the next map
        prefix.append("/.kmap." + std::to_string(m) + ".bin");
        
        phmap::flat_hash_map<uint64_t, VALUE> nextMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap.phmap_load(ar_in);
        
        unionSum(map[m], nextMap);
        
    }
    
    histogram(map[m]);
    
    freeContainer(map[m]);
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::unionSum(phmap::flat_hash_map<uint64_t, VALUE>& map1, phmap::flat_hash_map<uint64_t, VALUE>& map2) {
    
    for (auto pair : map2)
        map1[pair.first] += pair.second;
    
    return true;
    
}


template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::loadMap(std::string prefix, uint16_t m) { // loading prototype
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    map[m].phmap_load(ar_in);
    
    histogram(map[m]);
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    map[m].phmap_dump(ar_out);
    
    freeContainer(map[m]);
    
    return true;
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::report(INPUT& userInput) {
    
    const static phmap::flat_hash_map<std::string,int> string_to_case{
        {"stats",1},
        {"hist",2},
        {"kc",3}

    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    lg.verbose("Writing ouput: " + ext);
    
    // here we create a smart pointer to handle any kind of output stream
    std::unique_ptr<std::ostream> ostream;
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
            
        case 1: { // .stats
            
            std::ofstream ofs(userInput.outFile);
            
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            
            *ostream<<"\nTotal: "<<totKmers<<"\n";
            *ostream<<"Unique: "<<totKmersUnique<<"\n";
            *ostream<<"Distinct: "<<totKmersDistinct<<"\n";
            uint64_t missing = pow(4,k)-totKmersDistinct;
            *ostream<<"Missing: "<<missing<<"\n";
            
            ofs.close();
            
            break;
            
        }
        case 2: { // .hist
            
            std::ofstream ofs(userInput.outFile);
            
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            
            printHist(ostream);
            
            ofs.close();
            
            break;
            
        }
        case 3: { // .kc
            
            make_dir(userInput.outFile.c_str());
            
            std::ofstream ofs(userInput.outFile + "/.index");
            
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            
            *ostream<<+k<<"\n"<<mapCount<<std::endl;
            
            ofs.close();
            
            for(uint16_t m = 0; m<mapCount; ++m)
                threadPool.queueJob([=]{ return dumpMap(userInput.outFile, m); });
            
            jobWait(threadPool);
            
            break;
            
        }
        default: {
            
            ostream = std::make_unique<std::ostream>(std::cout.rdbuf());
            
            printHist(ostream);
            
            *ostream<<"\nTotal: "<<totKmers<<"\n";
            *ostream<<"Unique: "<<totKmersUnique<<"\n";
            *ostream<<"Distinct: "<<totKmersDistinct<<"\n";
            uint64_t missing = pow(4,k)-totKmersDistinct;
            *ostream<<"Missing: "<<missing<<"\n";
            
        }
            
    }
    
}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::traverseInReads(Sequences* readBatch) { // traverse the read

    hashSequences(readBatch);
    
    delete readBatch;
    
    return true;
    
}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::traverseInReads(std::string* readBatch) { // traverse the read

    hashSequences(readBatch);
    
    delete readBatch;
    
    return true;
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::appendReads(Sequences* readBatch) { // read a collection of reads
    
    threadPool.queueJob([=]{ return traverseInReads(readBatch); });
    
    std::unique_lock<std::mutex> lck(mtx);
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        
    }
    
}

template<class INPUT, typename VALUE, typename TYPE>
inline uint64_t Kmap<INPUT, VALUE, TYPE>::hash(uint8_t *kmer, bool *isFw) {
    
    uint64_t fw = 0, rv = 0;
    
    for(uint8_t c = 0; c<k; ++c)
        fw += *kmer++ * pows[c];
    
    --kmer;
    
    for(uint8_t c = 0; c<k; ++c)
        rv += (3-(*kmer--)) * pows[c];
    
    if (isFw != NULL)
        *isFw = fw < rv ? true : false;
    
    return fw < rv ? fw : rv;
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::consolidate() {
    
    for (unsigned int i = 0; i<buffers.size(); ++i) {
        
        unsigned int counter = 0;
        
        for(uint16_t m = 0; m<mapCount; ++m) {
            
            Buf<uint64_t>* thisBuf = &buffers[i][m];
            
            if (thisBuf->seq != NULL && mapsInUse[m] == false) {
                
                mapsInUse[m] = true;
                threadPool.queueJob([=]{ return countBuff(thisBuf, m); });
                
            }
            
            if(thisBuf->seq == NULL){
                
                ++counter;
                
                if (counter == mapCount) {
                    lg.verbose("Consolidating buffer " + std::to_string(i));
                    buffers.erase(buffers.begin() + i);
                }
                
            }

        }
        
    }

}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::finalize() {

    lg.verbose("Counting with " + std::to_string(mapCount) + " maps");
        
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return countBuffs(m); });

}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::countBuff(Buf<uint64_t>* thisBuf, uint16_t m) {

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    if (thisBuf->seq != NULL) {
        
        phmap::flat_hash_map<uint64_t, VALUE>* thisMap;
        
        thisMap = &map[m];
        
        uint64_t len = thisBuf->pos;
        
        for (uint64_t c = 0; c<len; ++c)
            ++(*thisMap)[thisBuf->seq[c]];
        
        delete[] thisBuf->seq;
        thisBuf->seq = NULL;
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    mapsInUse[m] = false;
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::countBuffs(uint16_t m) {

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    Buf<uint64_t>* thisBuf;
    
    phmap::flat_hash_map<uint64_t, VALUE>* thisMap;
    
    for(Buf<uint64_t>* buf : buffers) {
            
        thisBuf = &buf[m];
        
        if (thisBuf->seq != NULL) {
            
            thisMap = &map[m];
            
            uint64_t len = thisBuf->pos;
            
            for (uint64_t c = 0; c<len; ++c)
                ++(*thisMap)[thisBuf->seq[c]];
            
            delete[] thisBuf->seq;
            thisBuf->seq = NULL;
            
        }
        
    }
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::histogram(phmap::flat_hash_map<uint64_t, VALUE>& map) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0;
    
    phmap::flat_hash_map<uint64_t, VALUE> hist;
    
    for (auto pair : map) {
        
        if (pair.second == 1)
            ++kmersUnique;
        
        ++kmersDistinct;
        
        ++hist[pair.second];
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    totKmersUnique += kmersUnique;
    
    totKmersDistinct += kmersDistinct;
    
    for (auto pair : hist) {
        
        finalHistogram[pair.first] += pair.second;
        
        totKmers += pair.first * pair.second;
        
    }
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::printHist(std::unique_ptr<std::ostream>& ostream) {
    
    std::vector<std::pair<uint64_t, uint64_t>> table(finalHistogram.begin(), finalHistogram.end());
    std::sort(table.begin(), table.end());
    
    for (auto pair : table)
        *ostream<<pair.first<<"\t"<<pair.second<<"\n";

}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::hashSequences(Sequences* readBatch) {
    
    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    Buf<TYPE>* buf = new Buf<TYPE>[mapCount];
    
    for (Sequence* sequence : readBatch->sequences) {
        
        uint64_t len = sequence->sequence->size();
        
        if (len<k)
            continue;
        
        unsigned char* first = (unsigned char*)sequence->sequence->c_str();
        
        uint8_t* str = new uint8_t[len];
        uint64_t e = 0;
        
        for (uint64_t p = 0; p<len; ++p) {
            
            str[p] = ctoi[*(first+p)];
            
            if (str[p] > 3 || p+1==len){
                
                if (p+1==len && str[p] < 4) { // end of sequence, adjust indexes
                    ++e;
                    ++p;
                }
                
                if (e < k) { // beginning/end of a sequence or kmer too short, nothing to be done
                    e = 0;
                    continue;
                }
                
                uint64_t kcount = e-k+1;
                
                uint64_t key, i, newSize;
                Buf<TYPE>* b;
                TYPE* bufNew;
                
                for (uint64_t c = 0; c<kcount; ++c){
                    
                    key = hash(str+c+p-e);
                    i = key / moduloMap;
                    b = &buf[i];
                    
                    if (b->pos == b->size) {
                        
                        newSize = b->size * 2;
                        bufNew = new TYPE[newSize];
                        
                        memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));
                        
                        b->size = newSize;
                        delete[] b->seq;
                        b->seq = bufNew;
                        
                    }
                    
                    b->seq[b->pos++] = key;
                    
                }
                
                e = 0;
                
            }else{
                
                ++e;
                
            }
            
        }
        
        delete[] str;
        
//        threadLog.add("Processed sequence: " + sequence->header);
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
    logs.push_back(threadLog);
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::hashSequences(std::string* readBatch) {
    
    Log threadLog;
    
    Buf<TYPE>* buf = new Buf<TYPE>[mapCount];
        
    uint64_t len = readBatch->size();
    
    if (len<k)
        return;
    
    unsigned char* first = (unsigned char*)readBatch->c_str();
    
    uint8_t* str = new uint8_t[len];
    uint64_t e = 0;
    
    for (uint64_t p = 0; p<len; ++p) {
        
        str[p] = ctoi[*(first+p)];
        
        if (str[p] > 3 || p+1==len){
            
            if (p+1==len && str[p] < 4) { // end of sequence, adjust indexes
                ++e;
                ++p;
            }
            
            if (e < k) { // beginning/end of a sequence or kmer too short, nothing to be done
                e = 0;
                continue;
            }
            
            uint64_t kcount = e-k+1;
            
            uint64_t key, i, newSize;
            Buf<TYPE>* b;
            TYPE* bufNew;
            
            for (uint64_t c = 0; c<kcount; ++c){
                
                key = hash(str+c+p-e);
                i = key / moduloMap;
                b = &buf[i];
                
                if (b->pos == b->size) {
                    
                    newSize = b->size * 2;
                    bufNew = new TYPE[newSize];
                    
                    memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));
                    
                    b->size = newSize;
                    delete[] b->seq;
                    b->seq = bufNew;
                    
                }
                
                b->seq[b->pos++] = key;
                
            }
            
            e = 0;
            
        }else{
            
            ++e;
            
        }
        
    }
    
    delete[] str;
        
//        threadLog.add("Processed sequence: " + sequence->header);
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
    logs.push_back(threadLog);
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::hashSegments() {
    
    std::vector<InSegment*>* segments = inSequences.getInSegments();
    
    if (segments->size() == 0)
        return;
    
    Buf<TYPE>* buf = new Buf<TYPE>[mapCount];
    
    for (InSegment* segment : *segments) {
        
        uint64_t len = segment->getSegmentLen(), kcount = len-k+1;
        
        if (len<k)
            continue;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i)
            str[i] = ctoi[*(first+i)];
        
        uint64_t key, i, newSize;
        Buf<TYPE>* b;
        TYPE* bufNew;
        
        for (uint64_t c = 0; c<kcount; ++c){
            
            key = hash(str+c);
            
            i = key / moduloMap;
            
            b = &buf[i];
            
            if (b->pos == b->size) {
                
                newSize = b->size * 2;
                bufNew = new TYPE[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            b->seq[b->pos++] = key;
                        
        }
        
        delete[] str;
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::hist() {

    lg.verbose("Generate histogram");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return histogram(map[m]); });
    
    jobWait(threadPool);
    
}


#endif //KMER
