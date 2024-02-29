#ifndef KMER_H
#define KMER_H

#include <parallel-hashmap/phmap.h>
#include "parallel-hashmap/phmap_dump.h"

#include <fastx.h>

template<typename T>
inline void freeContainer(T& p_container) // this is a C++ trick to empty a container and release associated memory
{
    T empty;
    std::swap(p_container, empty); // swapping a container with an empty (NULL) container should release the associated memory
}

template<typename TYPE> // this is a generic buffer, VALUE is the type of the elements we wish to store in it. Usually each hashed kmer becomes part of a buffer specified by its hash value
struct Buf {
    uint64_t pos = 0, size; // pos keeps track of the position reached filling the buffer, initialized to contain up to size elements
    TYPE *seq = new TYPE[size]; // the actual container
    
    Buf() : size(pow(2,8)){
        alloc += size*sizeof(TYPE);
    }
    Buf(uint64_t size) : size(size){
        alloc += size*sizeof(TYPE);
    }
    
    uint64_t newPos(uint8_t bytes) {
        
        if (pos + bytes > size) {
            
            uint64_t newSize = size*2;
            alloc += newSize*sizeof(TYPE);
            TYPE* seqNew = new TYPE[newSize];
            
            memcpy(seqNew, seq, size*sizeof(TYPE));
            
            delete[] seq;
            freed += size*sizeof(TYPE);
            size = newSize;
            seq = seqNew;
            
        }
        
        return pos += bytes;
        
    }
};

template<class INPUT, typename VALUE, typename TYPE> // INPUT is a specialized userInput type depending on the tool, VALUE is the type of elements we wish to store in the maps, e.g. uint64_t kmer counts, TYPE is the type of inputs, e.g. kmers hashed to sequences
class Kmap {

protected: // they are protected, so that they can be further specialized by inheritance
    
    InSequences inSequences; // when we read a reference we can store it here
    
    uint32_t processedBuffers = 0; // useful to keep track of buffers as they are processed

    uint8_t k; // klen
    
    uint64_t totKmers = 0, totKmersUnique = 0, totKmersDistinct = 0; // summary statistics
    
    const uint16_t mapCount = 128; // number of maps to store the kmers, the longer the kmers, the higher number of maps to increase efficiency
    
    const uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount; // this value allows to assign any kmer to a map based on its hashed value
    
    uint64_t* pows = new uint64_t[k]; // storing precomputed values of each power significantly speeds up hashing

    std::vector<Buf<TYPE>*> buffersVec; // a vector for all buffers
    
    using parallelMap = phmap::parallel_flat_hash_map<uint64_t, VALUE,
                                              std::hash<uint64_t>,
                                              std::equal_to<uint64_t>,
                                              std::allocator<std::pair<const uint64_t, VALUE>>,
                                              8,
                                              phmap::NullMutex>;
    
    std::vector<parallelMap*> maps; // all hash maps where VALUES are stored
    
    std::vector<bool> mapsInUse = std::vector<bool>(mapCount, false); // useful with multithreading to ensure non-concomitant write access to maps
    
    phmap::flat_hash_map<uint64_t, uint64_t> finalHistogram; // the final kmer histogram
    
    const uint8_t ctoi[256] = { // this converts ACGT>0123 and any other character to 4 in time O(1)
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
    
    const uint8_t itoc[4] = {'A', 'C', 'G', 'T'}; // 0123>ACGT
    
public:
    
    Kmap(uint8_t k) : k(k) { // precomputes the powers of k
        
        for(uint8_t p = 0; p<k; ++p)
            pows[p] = (uint64_t) pow(4,p);
        
        for(uint16_t m = 0; m<mapCount; ++m)
            maps.push_back(new parallelMap);
        
    };
    
    ~Kmap(){ // always need to call the destructor and delete for any object called with new to avoid memory leaks
        
        for (parallelMap* map : maps)
            delete map;
        delete[] pows;
        
    }
    
    void load(INPUT& userInput);
    
    void kunion(INPUT& userInput);
    
    bool mergeMaps(std::vector<std::string> prefixes, uint16_t m);
    
    bool unionSum(parallelMap& map1, parallelMap& map2);
    
    bool traverseInReads(Sequences* readBatch);
    
    bool traverseInReads(std::string* readBatch);
    
    void appendSequence(Sequence* sequence, int hc_cutoff);
    
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
    
    bool histogram(parallelMap& map);
    
    void resizeBuff(Buf<uint64_t>* buff);
    
    void printHist(std::unique_ptr<std::ostream>& ostream);
    
    void report(INPUT& userInput);
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);

};

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::appendSequence(Sequence* sequence, int hc_cutoff) { // method to append a new sequence from a fasta
    
    threadPool.queueJob([=]{ return inSequences.traverseInSequence(sequence, hc_cutoff); }); // generic method to add a new job to the queue
    
    if(verbose_flag) {std::cerr<<"\n";};
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock); // every time a shared variable is edited we need to stop the threads
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::load(INPUT& userInput){ // concurrent loading of existing hashmaps
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return loadMap(userInput.iSeqFileArg, m); });
    
    jobWait(threadPool);
    
}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::loadMap(std::string prefix, uint16_t m) { // loads a specific maps
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    maps[m].phmap_load(ar_in);
    
    histogram(maps[m]);
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::kunion(INPUT& userInput){ // concurrent merging of the maps that store the same hashes
        
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return mergeMaps(userInput.iReadFileArg, m); });
    
    jobWait(threadPool);
    
}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::mergeMaps(std::vector<std::string> prefixes, uint16_t m) { // a single job merging maps with the same hashes
    
    std::string prefix = prefixes[0]; // loads the first map
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    maps[m].phmap_load(ar_in);
    
    unsigned int numFiles = prefixes.size();

    for (unsigned int i = 1; i < numFiles; i++) { // for each kmerdb loads the map and merges it
        
        std::string prefix = prefixes[i]; // loads the next map
        prefix.append("/.kmap." + std::to_string(m) + ".bin");
        
        parallelMap nextMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap.phmap_load(ar_in);
        
        unionSum(maps[m], nextMap); // unionSum operation between the existing map and the next map
        
    }
    
    histogram(maps[m]);
    
    freeContainer(maps[m]);
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::unionSum(parallelMap& map1, parallelMap& map2) {
    
    for (auto pair : map2) // for each element in map2, find it in map1 and increase its value
        map1[pair.first] += pair.second;
    
    return true;
    
}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m].phmap_dump(ar_out);
    
    freeContainer(maps[m]);
    
    return true;
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::report(INPUT& userInput) { // generates the output from the program
    
    const static phmap::parallel_flat_hash_map<std::string,int> string_to_case{ // different outputs available
        {"stats",1},
        {"hist",2},
        {"kc",3}

    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    lg.verbose("Writing ouput: " + ext);
    
    std::unique_ptr<std::ostream> ostream; // smart pointer to handle any kind of output stream
    
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
                threadPool.queueJob([=]{ return dumpMap(userInput.outFile, m); }); // writes map to file concurrently
            
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
bool Kmap<INPUT, VALUE, TYPE>::traverseInReads(Sequences* readBatch) { // specialized for Sequences objects

    hashSequences(readBatch);
    
    delete readBatch;
    
    return true;
    
}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::traverseInReads(std::string* readBatch) { // specialized for string objects

    hashSequences(readBatch);
    
    delete readBatch;
    
    return true;
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::appendReads(Sequences* readBatch) { // reads a collection of reads
    
    threadPool.queueJob([=]{ return traverseInReads(readBatch); });
    
}

template<class INPUT, typename VALUE, typename TYPE>
inline uint64_t Kmap<INPUT, VALUE, TYPE>::hash(uint8_t *kmer, bool *isFw) { // hashing function for kmers
    
    uint64_t fw = 0, rv = 0; // hashes for both forward and reverse complement sequence
    
    for(uint8_t c = 0; c<k; ++c) { // for each position up to klen
        fw += *kmer * pows[c]; // base * 2^N
        rv += (3-(*kmer++)) * pows[k-c-1]; // we walk the kmer backward to compute the rvcp
    }
    
    if (isFw != NULL)
        *isFw = fw < rv ? true : false; // we preserve the actual orientation for DBG applications
    
    return fw < rv ? fw : rv;
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::consolidate() { // to reduce memory footprint we consolidate the buffers as we go
    
    for (unsigned int i = 0; i<buffersVec.size(); ++i) { // for each buffer
        
        unsigned int counter = 0;
        
        for(uint16_t m = 0; m<mapCount; ++m) { // for each map
            
            Buf<uint64_t>* thisBuf = &buffersVec[i][m];
            
            if (thisBuf->seq != NULL && mapsInUse[m] == false) { // if the buffer was not counted and the associated map is not in use we process it
                
                mapsInUse[m] = true;
                threadPool.queueJob([=]{ return countBuff(thisBuf, m); });
                
            }
            
            if(thisBuf->seq == NULL){
                
                ++counter; // keeps track of the buffers that were processed so far
                
                if (counter == mapCount) {
                    lg.verbose("Jobs waiting/running: " + std::to_string(threadPool.queueSize()) + "/" + std::to_string(threadPool.running()) + " memory used/total: " + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
                    buffersVec.erase(buffersVec.begin() + i);
                }
                
            }

        }
        
    }

}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::finalize() { // ensure we count all residual buffers

    lg.verbose("Counting with " + std::to_string(mapCount) + " maps");
        
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return countBuffs(m); });

}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::countBuff(Buf<uint64_t>* thisBuf, uint16_t m) { // counts a single buffer

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    if (thisBuf->seq != NULL) { // sanity check that this buffer was not already processed
        
        parallelMap* thisMap = &maps[m]; // the map associated to this buffer
        
        uint64_t len = thisBuf->pos; // how many positions in the buffer have data
        
        for (uint64_t c = 0; c<len; ++c)
            ++(*thisMap)[thisBuf->seq[c]]; // writes to the map
        
        delete[] thisBuf->seq; // delete the buffer
        thisBuf->seq = NULL; // set its sequence to the null pointer in case its checked again
        
    }
    
    std::unique_lock<std::mutex> lck(mtx); // release the map
    mapsInUse[m] = false;
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
bool Kmap<INPUT, VALUE, TYPE>::countBuffs(uint16_t m) { // counts all residual buffers for a certain map as we finalize the kmerdb

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    Buf<uint64_t>* thisBuf;
    
    parallelMap* thisMap;
    
    for(Buf<uint64_t>* buf : buffersVec) {
            
        thisBuf = &buf[m];
        
        if (thisBuf->seq != NULL) {
            
            thisMap = &maps[m];
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
bool Kmap<INPUT, VALUE, TYPE>::histogram(parallelMap& map) { // extracts information from each map to build histogram
    
    uint64_t kmersUnique = 0, kmersDistinct = 0;
    
    parallelMap hist;
    
    for (auto pair : map) {
        
        if (pair.second == 1) // unique kmers
            ++kmersUnique;
        
        ++kmersDistinct; // distinct kmers
        ++hist[pair.second]; // increase the count of kmers with a certain frequency
        
    }
    
    std::unique_lock<std::mutex> lck(mtx); // updates the histogram, which is shared by all maps
    
    totKmersUnique += kmersUnique;
    
    totKmersDistinct += kmersDistinct;
    
    for (auto pair : hist) { // consolidate histograms from all maps into one
        
        finalHistogram[pair.first] += pair.second;
        totKmers += pair.first * pair.second; // number of kmers with a certain frequency times the frequency
        
    }
    
    return true;

}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::printHist(std::unique_ptr<std::ostream>& ostream) { // prints the histogram
    
    std::vector<std::pair<uint64_t, uint64_t>> table(finalHistogram.begin(), finalHistogram.end()); // converts the hashmap to a table
    std::sort(table.begin(), table.end());
    
    for (auto pair : table)
        *ostream<<pair.first<<"\t"<<pair.second<<"\n";

}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::hashSequences(Sequences* readBatch) { // hashes a batch of input reads
    
    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    Buf<TYPE>* buf = new Buf<TYPE>[mapCount]; // creates a temporary buffer for each map
    
    for (Sequence* sequence : readBatch->sequences) { // for each input read in this batch
        
        uint64_t len = sequence->sequence->size();
        
        if (len<k) // a read shorter than klen does not produce any kmer
            continue;
        
        unsigned char* first = (unsigned char*)sequence->sequence->c_str(); // useful point to the first character of the read
        
        uint8_t* str = new uint8_t[len]; // a new buffer to store the bases converted to numbers
        uint64_t e = 0; // contig length
        
        for (uint64_t p = 0; p<len; ++p) { // for each position (p) in the read
            
            str[p] = ctoi[*(first+p)]; // converts current base to number
            
            if (str[p] > 3 || p+1==len){ // we have found a gap or the sequence end
                
                if (p+1==len && str[p] < 4) { // end of sequence, adjust indexes
                    ++e;
                    ++p;
                }
                
                if (e < k) { // beginning/end of a sequence or kmer too short, nothing to be done
                    e = 0;
                    continue;
                }
                
                uint64_t kcount = e-k+1; // number of kmers in this contig
                
                uint64_t key, i, newSize;
                Buf<TYPE>* b; // a pointer to a buffer that is swapped based on the kmer hash value
                
                for (uint64_t c = 0; c<kcount; ++c){ // for each kmer in this contig
                    
                    key = hash(str+c+p-e); // hashes kmer
                    i = key / moduloMap; // finds the maps the kmer belongs to
                    b = &buf[i]; // assign the corresponding buffer
                    
                    if (b->pos == b->size) { // increases the size of the buffer if the end is reached
                        
                        newSize = b->size * 2;
                        TYPE* bufNew = new TYPE[newSize];
                        
                        memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));
                        
                        b->size = newSize;
                        delete[] b->seq;
                        b->seq = bufNew;
                        
                    }
                    
                    b->seq[b->pos++] = key; // write the hashed value in the buffer
                    
                }
                
                e = 0; // resets the contig length counter
                
            }else{ // increases contig length
                
                ++e;
                
            }
            
        }
        
        delete[] str; // deletes the numeric string
        
//        threadLog.add("Processed sequence: " + sequence->header);
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffersVec.push_back(buf); // stores all the new buffers just generated
    
    logs.push_back(threadLog);
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::hashSequences(std::string* readBatch) { // same as previous function, but when reads are stored as simple strings
    
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
    
    buffersVec.push_back(buf);
    
    logs.push_back(threadLog);
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::hashSegments() { // hashes (gapless) segments, no need to check for N bases
    
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
    
    buffersVec.push_back(buf);
    
}

template<class INPUT, typename VALUE, typename TYPE>
void Kmap<INPUT, VALUE, TYPE>::hist() {

    lg.verbose("Generate histogram");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return histogram(maps[m]); });
    
    jobWait(threadPool);
    
}


#endif //KMER
