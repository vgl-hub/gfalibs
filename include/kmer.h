#ifndef KMER_H
#define KMER_H

#include <future>
#include <parallel-hashmap/phmap.h>
#include "parallel-hashmap/phmap_dump.h"

#include <fastx.h>

#define LARGEST 4294967295 // 2^32-1

template<typename T>
inline void freeContainer(T& p_container) // this is a C++ trick to empty a container and release associated memory
{
    T empty;
    std::swap(p_container, empty); // swapping a container with an empty (NULL) container should release associated memory
}

class Key {
    uint64_t kmer;
public:
    
    Key(){}
    Key(uint64_t kmer) : kmer(kmer){}
    
    void assignKey(uint64_t kmer) {
        this->kmer = kmer;
    }
    uint64_t getKmer() const {
        return kmer;
    }
};

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2> // DERIVED implements the CRTP technique, INPUT is a specialized userInput type depending on the tool, KEY is the key type for the hashtable, TYPE1 is the low frequency type of elements we wish to store in the maps, e.g. uint8_t kmer counts, TYPE2 is the high frequency type of elements we wish to store in the maps, e.g. uint32_t kmer counts
class Kmap {

protected: // they are protected, so that they can be further specialized by inheritance
    
    UserInput &userInput;
    InSequences inSequences; // when we read a reference we can store it here
    
    uint8_t k; // kPrefixLen
    uint64_t tot = 0, totUnique = 0, totDistinct = 0; // summary statistics
    std::atomic<bool> readingDone{false};
    std::vector<std::thread> threads;
    std::vector<std::future<bool>> futures;
    std::string DBextension;
    
    const uint16_t mapCount = 127; // number of maps to store the kmers, the longer the kmers, the higher number of maps to increase efficiency
    
    const uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount; // this value allows to assign any kmer to a map based on its hashed value
    
    uint64_t* pows = new uint64_t[k]; // storing precomputed values of each power significantly speeds up hashing

    std::vector<Buf<uint64_t>*> buffersVec; // a vector for all kmer buffers
    std::vector<Buf<uint8_t>*> strVec; // sequence vector
    uint64_t strVecLen = 0;
    
    struct KeyHasher {

        uint64_t* pows = new uint64_t[kPrefixLen];
        
        KeyHasher() {
            for(uint8_t p = 0; p<kPrefixLen; ++p) // precomputes the powers of k
                pows[p] = (uint64_t) pow(4,p);
        }
        
        std::size_t operator()(const Key& key) const {
            
            uint8_t *kmerPtr = seqBuf->seq+key.getKmer();
            uint64_t fw = 0, rv = 0; // hashes for both forward and reverse complement sequence
            
            for(uint8_t c = 0; c<kPrefixLen; ++c) { // for each position up to kPrefixLen
                fw += *(kmerPtr+c) * pows[c]; // base * 2^N
                rv += (3-(*(kmerPtr+kLen-1-c))) * pows[c]; // we walk the kmer backward to compute the rvcp
            }
            
            // return fw < rv ? 0 : 0; // even if they end up in the same bucket it's fine!
            return fw < rv ? fw : rv;
        }
    };
    
    struct KeyEqualTo {
        
        constexpr bool operator()(const Key& key1, const Key& key2) const {
            
            uint8_t *lhs = seqBuf->seq+key1.getKmer(), *rhs = seqBuf2->seq+key2.getKmer();
            
            for(uint32_t i = 0; i<kLen; ++i) { // check fw
                if(lhs[i] ^ rhs[i])
                    break;
                if (i == kLen-1)
                    return true;
            }
            for(uint32_t i = 0; i<kLen; ++i) { // if fw fails, check rv
                if(lhs[i] ^ (3-rhs[kLen-i-1]))
                    return false;
            }
            return true;
        }
    };
    
    using ParallelMap = phmap::parallel_flat_hash_map<KEY, TYPE1,
                                              KeyHasher,
                                              KeyEqualTo,
                                              std::allocator<std::pair<const KEY, TYPE1>>,
                                              8,
                                              phmap::NullMutex>;
    
    using ParallelMap32 = phmap::parallel_flat_hash_map<KEY, TYPE2,
                                              KeyHasher,
                                              KeyEqualTo,
                                              std::allocator<std::pair<const KEY, TYPE2>>,
                                              8,
                                              phmap::NullMutex>;
    
    std::vector<ParallelMap*> maps; // all hash maps where TYPE1S are stored
    std::vector<ParallelMap32*> maps32;
    
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
    
    std::chrono::high_resolution_clock::time_point past;
    std::queue<std::string*> readBatches;
    std::queue<Sequences*> sequenceBatches;
    
    std::mutex readMtx, hashMtx;
    
public:
    
    Kmap(UserInput& userInput) : userInput(userInput), k{kPrefixLen} {
        
        DBextension = "kc";
        
        seqBuf = new Buf<uint8_t>;
        seqBuf2 = seqBuf;
        
        for(uint8_t p = 0; p<k; ++p) // precomputes the powers of k
            pows[p] = (uint64_t) pow(4,p);
        
        for(uint16_t m = 0; m<mapCount; ++m)
            maps.push_back(new ParallelMap);
        
        for(uint16_t m = 0; m<mapCount; ++m)
            maps32.push_back(new ParallelMap32);
        
        if (userInput.kmerDB.size() == 0) { // if we are not reading an existing db
            lg.verbose("Deleting any tmp file");
            for(uint16_t m = 0; m<mapCount; ++m) {// remove tmp buffers and maps if any
                threadPool.queueJob([=]{ return remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str()); });
                threadPool.queueJob([=]{ return remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str()); });
                uint8_t fileNum = 0;
                while (fileExists(userInput.prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum++) +  ".tmp.bin"))
                    threadPool.queueJob([=]{ return remove((userInput.prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin").c_str()); });
                remove((userInput.prefix + "/.index").c_str());
                remove((userInput.prefix + "/.seq.bin").c_str());

            }
            jobWait(threadPool);
            initHashing(); // start parallel hashing
        }
        
    };
    
    ~Kmap(){ // always need to call the destructor and delete for any object called with new to avoid memory leaks
        for (ParallelMap* map : maps)
            delete map;
        for (ParallelMap32* map : maps32)
            delete map;
        delete[] pows;
        delete seqBuf;
    }
    
    uint64_t mapSize(ParallelMap& m);
    
    void joinThreads();
    
    bool memoryOk();
    
    bool memoryOk(int64_t delta);
    
    void initHashing();
    
    bool dumpBuffers();
    
    void buffersToMaps();
    
    bool processBuffers(uint16_t m);
    
    void consolidateTmpMaps();
    
    bool dumpTmpMap(std::string prefix, uint16_t m);
    
    bool deleteMap(uint16_t m);
    
    void dumpHighCopyKmers();
    
    bool mergeTmpMaps(uint16_t m);
    
    bool reloadMap32(uint16_t m);
    
    void status();
    
    void kunion();
    
    bool mergeSubMaps(ParallelMap* map1, ParallelMap* map2, uint8_t subMapIndex, uint16_t m);
    
    bool unionSum(ParallelMap* map1, ParallelMap* map2, uint16_t m);
    
    bool traverseInReads(std::string* readBatch);
    
    bool traverseInReads(Sequences* readBatch);
    
    inline uint64_t hash(uint8_t* string, bool* isFw = NULL);
    
    inline std::string reverseHash(uint64_t hash);
    
    bool hashSequences();
    
    void consolidate();
    
    void finalize();
    
    void stats();
    
    void DBstats();
    
    bool summary(uint16_t m);
    
    void printHist(std::unique_ptr<std::ostream>& ostream);
    
    void report();
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);
    
    bool loadHighCopyKmers();
    
    std::array<uint16_t, 2> computeMapRange(std::array<uint16_t, 2> mapRange);
    
    void loadMapRange(std::array<uint16_t, 2> mapRange);
    
    void deleteMapRange(std::array<uint16_t, 2> mapRange);
    
    void cleanup();
    
    bool mergeMaps(uint16_t m);

};


template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::memoryOk() {
    
    return get_mem_inuse(3) < maxMem;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::memoryOk(int64_t delta) {
    
    return get_mem_inuse(3) + convert_memory(delta, 3) < maxMem;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
uint64_t Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::mapSize(ParallelMap& m) {
    
   return m.capacity() * (sizeof(typename ParallelMap::value_type) + 1) + sizeof(ParallelMap);
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::joinThreads() {
    
    uint8_t threadsDone = 0;
    bool done = false;
    
    while (!done) {
        for (uint8_t i = 0; i < futures.size(); ++i) {
            if (futures[i].wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
                if (threads[i].joinable()) {
                    threads[i].join();
                    ++threadsDone;
                }
            }else{
                status();
            }
        }
        if (threadsDone == futures.size())
            done = true;
    }
    futures.clear();
    threads.clear();
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::initHashing(){
    
    std::packaged_task<bool()> task([this] { return dumpBuffers(); });
    futures.push_back(task.get_future());
    threads.push_back(std::thread(std::move(task)));
    
    int16_t threadN = threadPool.totalThreads() - 1; // substract the writing thread
    
    if (threadN == 0)
        threadN = 1;
    
    for (uint8_t t = 0; t < threadN; t++) {
        std::packaged_task<bool()> task([this] { return static_cast<DERIVED*>(this)->hashSequences(); });
        futures.push_back(task.get_future());
        threads.push_back(std::thread(std::move(task)));
    }
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::dumpBuffers() {
    
    bool hashing = true;
    std::vector<Buf<uint64_t>*> buffersVecCpy;
    std::ofstream bufFiles[mapCount];
    
    for (uint16_t b = 0; b<mapCount; ++b) // we open all files at once
        bufFiles[b] = std::ofstream(userInput.prefix + "/.buf." + std::to_string(b) + ".bin", std::fstream::app | std::ios::out | std::ios::binary);
    
    while (hashing) {
        
        if (readingDone) { // if we have finished reading the input
            
            uint8_t hashingDone = 0;
        
            for (uint8_t i = 1; i < futures.size(); ++i) { // we check how many hashing threads are still running
                if (futures[i].wait_for(std::chrono::milliseconds(0)) == std::future_status::ready)
                    ++hashingDone;
            }
            
            if (hashingDone == futures.size() - 1) // if all hashing threads are done we exit the loop after one more iteration
                hashing = false;
        }
        
        {
            std::unique_lock<std::mutex> lck(hashMtx); // we safely collect the new buffers
            buffersVecCpy = buffersVec;
            buffersVec.clear();
        }
        
        for (Buf<uint64_t>* buffers : buffersVecCpy) { // for each array of buffers
            
            for (uint16_t b = 0; b<mapCount; ++b) { // for each buffer file
                Buf<uint64_t>* buffer = &buffers[b];
                bufFiles[b].write(reinterpret_cast<const char *>(&buffer->pos), sizeof(uint64_t));
                bufFiles[b].write(reinterpret_cast<const char *>(buffer->seq), sizeof(uint64_t) * buffer->pos);
                delete[] buffer->seq;
                buffer->seq = NULL;
                freed += buffer->size * sizeof(uint64_t);
            }
            delete[] buffers;
        }
    }
    
    for (uint16_t b = 0; b<mapCount; ++b) // we close all files
        bufFiles[b].close();
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::buffersToMaps() {
    
    std::vector<std::function<bool()>> jobs;
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of buf files
        fileSizes.push_back(fileSize(userInput.prefix + "/.buf." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        jobs.push_back([this, i] { return static_cast<DERIVED*>(this)->processBuffers(i); });
        
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
    
    consolidateTmpMaps();
    dumpHighCopyKmers();

}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::processBuffers(uint16_t m) {
    
    uint64_t pos = 0;
    Buf<uint64_t> *buf;
    
    std::string fl = userInput.prefix + "/.buf." + std::to_string(m) + ".bin";
    std::ifstream bufFile(fl, std::ios::in | std::ios::binary);
    
    while(bufFile && !(bufFile.peek() == EOF)) {
        
        {
            std::unique_lock<std::mutex> lck(mtx);
            std::condition_variable &mutexCondition = threadPool.getMutexCondition();
            mutexCondition.wait(lck, [] {
                return !freeMemory;
            });
        }
        
        ParallelMap& map = *maps[m]; // the map associated to this buffer
        ParallelMap32& map32 = *maps32[m];
        uint64_t map_size = mapSize(map);
        
        bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
        
        buf = new Buf<uint64_t>(pos);
        
        buf->pos = pos;
        buf->size = pos;
        
        bufFile.read(reinterpret_cast<char *>(buf->seq), sizeof(uint64_t) * buf->pos);
        
        for (uint64_t c = 0; c<pos; ++c) {
            
            Key key(buf->seq[c]);

            uint8_t &count = map[key];
            bool overflow = (count >= 254 ? true : false);
            
            if (!overflow)
                ++count; // increase kmer coverage
            
            if (overflow) {
                
                uint32_t &count32 = map32[key];
                
                if (count32 == 0) { // first time we add the kmer
                    count32 = count;
                    count = 255; // invalidates int8 kmer
                }

                if (count32 < LARGEST)
                    ++count32; // increase kmer coverage
            }
        }
        
        freed += buf->size * sizeof(uint8_t);
        delete buf;
        alloc += mapSize(*maps[m]) - map_size;
        
        if (freeMemory || !bufFile || bufFile.peek() == EOF) {
            dumpTmpMap(userInput.prefix, m);
            reloadMap32(m);
        }
    }

    bufFile.close();
    remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str());
    
    return true;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::reloadMap32(uint16_t m) {
    
    ParallelMap& map = *maps[m]; // the map associated to this buffer
    ParallelMap32& map32 = *maps32[m];
    
    for (auto pair : map32) {
        
        TYPE1 count = 255;
        auto newPair = std::make_pair(pair.first, count);
        map.insert(newPair);
    }
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::mergeTmpMaps(uint16_t m) { // a single job merging maps with the same hashes
    
    std::string prefix = userInput.prefix; // loads the first map
    std::string firstFile = prefix + "/.map." + std::to_string(m) + ".0.tmp.bin";
    
    if (!fileExists(prefix + "/.map." + std::to_string(m) + ".1.tmp.bin")) {
        
        std::rename(firstFile.c_str(), (prefix + "/.map." + std::to_string(m) + ".bin").c_str());
        return true;
        
    }
    
    uint8_t fileNum = 0;
    
    while (fileExists(prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) + ".tmp.bin")) { // for additional map loads the map and merges it
        
        std::string nextFile = prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum++) + ".tmp.bin"; // loads the next map
        ParallelMap* nextMap = new ParallelMap;
        phmap::BinaryInputArchive ar_in(nextFile.c_str());
        nextMap->phmap_load(ar_in);
        uint64_t map_size1 = mapSize(*nextMap);
        alloc += map_size1;
        
        uint64_t map_size2 = mapSize(*maps[m]);
        unionSum(nextMap, maps[m], m); // unionSum operation between the existing map and the next map
        
        alloc += mapSize(*maps[m]) - map_size2;
        remove(nextFile.c_str());
        delete nextMap;
        freed += map_size1;
        
    }
    dumpMap(userInput.prefix, m);
    
    return true;

}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::consolidateTmpMaps(){ // concurrent merging of the maps that store the same hashes
    
    lg.verbose("Consolidating temporary maps");
    
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
        fileSizes.push_back(fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".0.tmp.bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        mergeTmpMaps(i);
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::dumpTmpMap(std::string prefix, uint16_t m) {
    
    uint8_t fileNum = 0;
    
    while (fileExists(prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin"))
        ++fileNum;
        
    prefix.append("/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m]->phmap_dump(ar_out);
    
    deleteMap(m);
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::deleteMap(uint16_t m) {
    
    uint64_t map_size = mapSize(*maps[m]);
    delete maps[m];
    freed += map_size;
    
    maps[m] = new ParallelMap;
    alloc += mapSize(*maps[m]);
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::dumpHighCopyKmers() {
    
    ParallelMap32 map32Total;
    
    for (uint16_t m = 0; m<mapCount; ++m) {
        for (auto pair : *maps32[m])
            map32Total.insert(pair);
        delete maps32[m];
        maps32[m] = new ParallelMap32;
    }

    phmap::BinaryOutputArchive ar_out((userInput.prefix + "/.map.hc.bin").c_str());
    map32Total.phmap_dump(ar_out);
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::status() {
    
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - past;
    
    if (elapsed.count() > 0.1) {
        lg.verbose("Read batches: " + std::to_string(readBatches.size() + sequenceBatches.size()) + ". Hash buffers: " + std::to_string(buffersVec.size()) + ". Memory in use/allocated/total: " + std::to_string(get_mem_inuse(3)) + "/" + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
    
        past = std::chrono::high_resolution_clock::now();
    }
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::loadMap(std::string prefix, uint16_t m) { // loads a specific map
    
    prefix.append("/.map." + std::to_string(m) + ".bin");
    phmap::BinaryInputArchive ar_in(prefix.c_str());
//    allocMemory(fileSize(prefix));
    maps[m]->phmap_load(ar_in);
    alloc += mapSize(*maps[m]);
    
    return true;

}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::loadHighCopyKmers() {
    
    ParallelMap32 map32Total;
    phmap::BinaryInputArchive ar_in((userInput.prefix + "/.map.hc.bin").c_str());
    map32Total.phmap_load(ar_in);
    
    for (auto pair : map32Total) {
//        uint8_t i = hash(seqBuf->seq+pair.first.getKmer()) % mapCount;
        std::cout<<+pair.first.getKmer()<<std::endl;
//        maps32[i]->insert(pair);
    }
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::kunion(){ // concurrent merging of the maps that store the same hashes
    
    ParallelMap32 map32Total; // first merge high-copy kmers
    
    for (unsigned int i = 0; i < userInput.kmerDB.size(); ++i) { // for each kmerdb loads the map and merges it
        
        std::string prefix = userInput.kmerDB[i]; // loads the next map
        prefix.append("/.map.hc.bin");
        
        ParallelMap32 nextMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap.phmap_load(ar_in);
        
        for (auto pair : nextMap) {
            
            TYPE2& count32 = map32Total[pair.first];
            
            if (LARGEST - count32 >= pair.second)
                count32 += pair.second; // increase kmer coverage
            else
                count32 = LARGEST;
        }
    }
    
    for (auto pair : map32Total) {
        uint64_t i = hash(seqBuf->seq+pair.first.getKmer()) % mapCount;
        maps32[i]->insert(pair);
    }
    
    std::vector<std::function<bool()>> jobs;
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
        fileSizes.push_back(fileSize(userInput.kmerDB[0] + "/.map." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        static_cast<DERIVED*>(this)->mergeMaps(i);
    
    dumpHighCopyKmers();
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::mergeMaps(uint16_t m) { // a single job merging maps with the same hashes
    
    std::string prefix = userInput.kmerDB[0]; // loads the first map
    prefix.append("/.map." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    maps[m]->phmap_load(ar_in);
    
    uint64_t pos;
    std::ifstream bufFile = std::ifstream(userInput.kmerDB[0]+ "/.seq.bin", std::ios::in | std::ios::binary);
    bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
    seqBuf = new Buf<uint8_t>(pos);
    seqBuf->pos = pos;
    seqBuf->size = pos;
    bufFile.read(reinterpret_cast<char *>(seqBuf->seq), sizeof(uint8_t) * seqBuf->pos);

    for (unsigned int i = 1; i < userInput.kmerDB.size(); ++i) { // for each kmerdb loads the map and merges it
        
        bufFile = std::ifstream(userInput.kmerDB[0] + "/.seq.bin", std::ios::in | std::ios::binary);
        bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
        seqBuf2 = new Buf<uint8_t>(pos);
        seqBuf2->pos = pos;
        seqBuf2->size = pos;
        bufFile.read(reinterpret_cast<char *>(seqBuf2->seq), sizeof(uint8_t) * seqBuf2->pos);
        
        std::string prefix = userInput.kmerDB[i]; // loads the next map
        prefix.append("/.map." + std::to_string(m) + ".bin");
        
        ParallelMap* nextMap = new ParallelMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap->phmap_load(ar_in);
        
        unionSum(nextMap, maps[m], m); // unionSum operation between the existing map and the next map
        delete nextMap;
        delete seqBuf2;
    }
    
    seqBuf2 = seqBuf;
    
    dumpMap(userInput.prefix, m);
    deleteMap(m);
    
    static_cast<DERIVED*>(this)->summary(m);
    
    return true;

}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.map." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m]->phmap_dump(ar_out);
    
    deleteMap(m);
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::report() { // generates the output from the program
    
    const static phmap::parallel_flat_hash_map<std::string,int> string_to_case{ // different outputs available
        {"hist",1},
        {"kc",2}

    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    static_cast<DERIVED*>(this)->stats();
    
    lg.verbose("Writing ouput: " + ext);
    
    std::unique_ptr<std::ostream> ostream; // smart pointer to handle any kind of output stream
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
            
        case 1: { // .hist
            std::ofstream ofs(userInput.outFile);
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            printHist(ostream);
            ofs.close();
            break;
            
        }
        case 2: { // .kc
            std::ofstream ofs(userInput.outFile + "/.index"); // adding index
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            *ostream<<+kLen<<"\n"<<mapCount<<std::endl;
            ofs.close();
            break;
        }
        default: {}
    }
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::traverseInReads(std::string* readBatch) { // specialized for string objects
    
    while(freeMemory) {status();}
    
    {
        std::lock_guard<std::mutex> lck(readMtx);
        readBatches.push(readBatch);
    }
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::traverseInReads(Sequences* sequenceBatch) { // specialized for sequence objects
    
    while(freeMemory) {status();}
    
    {
        std::lock_guard<std::mutex> lck(readMtx);
        sequenceBatches.push(sequenceBatch);
    }
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
inline uint64_t Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::hash(uint8_t *kmerPtr, bool *isFw) { // hashing function for kmers
    
    uint64_t fw = 0, rv = 0; // hashes for both forward and reverse complement sequence
    
    for(uint8_t c = 0; c<k; ++c) { // for each position up to kPrefixLen
        fw += *(kmerPtr+c) * pows[c]; // base * 2^N
        rv += (3-(*(kmerPtr+kLen-1-c))) * pows[c]; // we walk the kmer backward to compute the rvcp
    }
    
    if (isFw != NULL)
        *isFw = fw < rv ? true : false; // we preserve the actual orientation for DBG applications
    
    return fw < rv ? fw : rv;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
inline std::string Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::reverseHash(uint64_t hash) { // hashing function for kmers
    
    std::string seq(k, 'A');
    
    for(uint8_t c = k; c > 0; --c) { // for each position up to klen
        uint8_t i = c-1; // to prevent overflow
        seq[i] = itoc[hash / pows[i]]; // base * 2^N
        hash = hash % pows[i]; // we walk the kmer backward to compute the rvcp
    }
    
//    if (hash != 0) {
//        std::cout<<"reashing error!"<<std::endl;
//        exit(EXIT_FAILURE);
//    }
    
    return seq;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::consolidate() { // to reduce memory footprint we consolidate the buffers as we go
    status();
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::finalize() { // ensure we count all residual buffers
    if (userInput.kmerDB.size() == 0) {
        readingDone = true;
        joinThreads();

        lg.verbose("Writing compressed sequences to disk");
        std::ofstream bufFileW = std::ofstream(userInput.prefix + "/.seq.bin", std::fstream::app | std::ios::out | std::ios::binary);
        bufFileW.write(reinterpret_cast<const char *>(&strVecLen), sizeof(uint64_t));
        for (Buf<uint8_t>* str : strVec) { // for each sequence
            bufFileW.write(reinterpret_cast<const char *>(str->seq), sizeof(uint8_t) * str->pos);
            delete str;
        }
        bufFileW.close();
        
        lg.verbose("Reading compressed sequences to memory");
        std::ifstream bufFileR = std::ifstream(userInput.prefix + "/.seq.bin", std::ios::in | std::ios::binary);
        delete seqBuf;
        bufFileR.read(reinterpret_cast<char *>(&strVecLen), sizeof(uint64_t));
        seqBuf = new Buf<uint8_t>(strVecLen);
        bufFileR.read(reinterpret_cast<char *>(seqBuf->seq), sizeof(uint8_t) * strVecLen);
        seqBuf->pos = strVecLen;
        bufFileR.close();
        seqBuf2 = seqBuf;
        
        lg.verbose("Converting buffers to maps");
        static_cast<DERIVED*>(this)->buffersToMaps();
        
        loadHighCopyKmers();
    }
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::stats() {
    
    lg.verbose("Computing summary statistics");
    
    std::vector<std::function<bool()>> jobs;
    std::array<uint16_t, 2> mapRange = {0,0};
    
    while (mapRange[1] < mapCount) {
        
        mapRange = computeMapRange(mapRange);
        loadMapRange(mapRange);
        
        for (uint32_t i = mapRange[0]; i < mapRange[1]; ++i)
            jobs.push_back([this, i] { return static_cast<DERIVED*>(this)->summary(i); });
        
        threadPool.queueJobs(jobs);
        jobWait(threadPool);
        jobs.clear();
        deleteMapRange(mapRange);
    }
    static_cast<DERIVED*>(this)->DBstats();
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::DBstats() {
    
    uint64_t missing = pow(4,kLen)-totDistinct;
    
    std::cout<<"DB Summary statistics:\n"
             <<"Total kmers: "<<tot<<"\n"
             <<"Unique kmers: "<<totUnique<<"\n"
             <<"Distinct kmers: "<<totDistinct<<"\n"
             <<"Missing kmers: "<<missing<<"\n";
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::summary(uint16_t m) {
    
    uint64_t unique = 0, distinct = 0;
    phmap::parallel_flat_hash_map<uint64_t, uint64_t> hist;
    
    for (auto pair : *maps[m]) {
        
        if (pair.second == 255) // check the large table
            continue;
        
        if (pair.second == 1)
            ++unique;
        
        ++distinct;
        ++hist[pair.second];
    }
    
    for (auto pair : *maps32[m]) {
        
        ++distinct;
        ++hist[pair.second];
    }
 
    std::lock_guard<std::mutex> lck(mtx);
    totUnique += unique;
    totDistinct += distinct;
    
    for (auto pair : hist) {
        
        finalHistogram[pair.first] += pair.second;
        tot += pair.first * pair.second;
    }
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::printHist(std::unique_ptr<std::ostream>& ostream) { // prints the histogram
    
    std::vector<std::pair<uint64_t, uint64_t>> table(finalHistogram.begin(), finalHistogram.end()); // converts the hashmap to a table
    std::sort(table.begin(), table.end());
    
    for (auto pair : table)
        *ostream<<pair.first<<"\t"<<pair.second<<"\n";

}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::hashSequences() {
    //   Log threadLog;
    std::string *readBatch;
    uint64_t len = 0, currentPos = 0;
    Buf<uint8_t> *str;
    
    while (true) {
            
        {
            std::unique_lock<std::mutex> lck(readMtx);
            
            if (readingDone && readBatches.size() == 0)
                return true;

            if (readBatches.size() == 0)
                continue;
            
            std::condition_variable &mutexCondition = threadPool.getMutexCondition();
            mutexCondition.wait(lck, [] {
                return !freeMemory;
            });
            
            readBatch = readBatches.front();
            readBatches.pop();
            len = readBatch->size();
            
            if (len<kLen) {
                delete readBatch;
                freed += len * sizeof(char);
                continue;
            }
            str = new Buf<uint8_t>(len);
            str->pos = len;
            strVec.push_back(str);
            currentPos = strVecLen;
            strVecLen += len;
        }

        Buf<uint64_t> *buffers = new Buf<uint64_t>[mapCount];
        const unsigned char *first = (unsigned char*) readBatch->c_str();
        uint32_t e = 0;
        uint64_t key, pos = 0, kcount = len-kLen+1;
        bool isFw = false;
        Buf<uint64_t>* buffer;
        
        for (uint64_t p = 0; p<kcount; ++p) {
            
            for (uint32_t c = e; c<kLen; ++c) { // generate k bases if e=0 or the next if e=k-1
                
                str->seq[p+c] = ctoi[*(first+p+c)]; // convert the next base
                if (str->seq[p+c] > 3) { // if non-canonical base is found
                    p = p+c; // move position
                    e = 0; // reset base counter
                    break;
                }
                e = kLen-1;
            }
            
            if (e == 0) // not enough bases for a kmer
                continue;

            key = hash(&str->seq[p], &isFw);

            buffer = &buffers[key % mapCount];
            pos = buffer->newPos(1);
            
            buffer->seq[pos-1] = currentPos+p;
        }
        
        delete readBatch;
        //    threadLog.add("Processed sequence: " + sequence->header);
        //    std::lock_guard<std::mutex> lck(mtx);
        //    logs.push_back(threadLog);
        std::lock_guard<std::mutex> lck(hashMtx);
        freed += len * sizeof(char);
        buffersVec.push_back(buffers);
    }
    return true;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::mergeSubMaps(ParallelMap* map1, ParallelMap* map2, uint8_t subMapIndex, uint16_t m) {
    
    auto& inner = map1->get_inner(subMapIndex);   // to retrieve the submap at given index
    auto& submap1 = inner.set_;        // can be a set or a map, depending on the type of map1
    auto& inner2 = map2->get_inner(subMapIndex);
    auto& submap2 = inner2.set_;
    ParallelMap32& map32 = *maps32[m];
    
    for (auto pair : submap1) { // for each element in map1, find it in map2 and increase its value
        
        bool overflow = false;
        
        if (pair.second == 255) // already added to int32 map
            continue;
        
        auto got = map32.find(pair.first); // check if this is already a high-copy kmer
        if (got != map32.end()) {
            overflow = true;
        }else{
            
            auto got = submap2.find(pair.first); // insert or find this kmer in the hash table
            if (got == submap2.end()) {
                submap2.insert(pair);
            }else{
                
                TYPE1& count = got->second;
                    
                if (255 - count <= pair.second)
                    overflow = true;
                
                if (!overflow)
                    count += pair.second; // increase kmer coverage
            }
        }
        
        if (overflow) {
            
            TYPE2& count32 = map32[pair.first];
            
            if (count32 == 0) { // first time we add the kmer
                auto got = submap2.find(pair.first);
                TYPE1& count = got->second;
                count32 = count;
                count = 255; // invalidates int8 kmer
            }
            
            if (LARGEST - count32 >= pair.second)
                count32 += pair.second; // increase kmer coverage
            else
                count32 = LARGEST;
        }
    }
    return true;
}


template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::unionSum(ParallelMap* map1, ParallelMap* map2, uint16_t m) {
    
    std::vector<std::function<bool()>> jobs;
    
    if (map1->subcnt() != map2->subcnt()) {
        fprintf(stderr, "Maps don't have the same numbers of submaps (%zu != %zu). Terminating.\n", map1->subcnt(), map2->subcnt());
        exit(EXIT_FAILURE);
    }
    
    for(std::size_t subMapIndex = 0; subMapIndex < map1->subcnt(); ++subMapIndex)
        jobs.push_back([this, map1, map2, subMapIndex, m] { return static_cast<DERIVED*>(this)->mergeSubMaps(map1, map2, subMapIndex, m); });
    
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
    
    return true;
    
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
std::array<uint16_t, 2> Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::computeMapRange(std::array<uint16_t, 2> mapRange) {
    
    uint64_t max = 0;
    mapRange[0] = mapRange[1];
    
    for (uint16_t m = mapRange[0]; m<mapCount; ++m) {
        
        max += fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".bin");
        if(!memoryOk(max))
            break;
        mapRange[1] = m + 1;
        
    }
    return mapRange;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::loadMapRange(std::array<uint16_t, 2> mapRange) {
    
    std::vector<std::function<bool()>> jobs;
    
    for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
        jobs.push_back([this, m] { return loadMap(userInput.prefix, m); });
    
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::deleteMapRange(std::array<uint16_t, 2> mapRange) {
    
    for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
        deleteMap(m);
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::cleanup() {
    
    if(!(userInput.kmerDB.size() == 1) && userInput.outFile.find("." + DBextension) == std::string::npos) {
        
        lg.verbose("Deleting tmp files");
        
        std::vector<std::function<bool()>> jobs;
        
        for(uint16_t m = 0; m<mapCount; ++m) // remove tmp files
            jobs.push_back([this, m] { return remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str()); });
        
        threadPool.queueJobs(jobs);
        
        remove((userInput.prefix + "/.map.hc.bin").c_str());
        remove((userInput.prefix + "/.seq.bin").c_str());
        
        if (userInput.prefix != ".")
            rm_dir(userInput.prefix.c_str());
        
    }
    
    jobWait(threadPool);
    
}

#endif //KMER
