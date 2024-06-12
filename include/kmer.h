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
    std::swap(p_container, empty); // swapping a container with an empty (NULL) container should release the associated memory
}

template<typename TYPE> // this is a generic buffer, TYPE is the type of the elements we wish to store in it. Usually each hashed kmer becomes part of a buffer specified by its hash value
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

template<class INPUT, typename TYPE1, typename TYPE2> // INPUT is a specialized userInput type depending on the tool, TYPE21 is the low frequency type of elements we wish to store in the maps, e.g. uint8_t kmer counts, TYPE22 is the high frequency type of elements we wish to store in the maps, e.g. uint32_t kmer counts
class Kmap {

protected: // they are protected, so that they can be further specialized by inheritance
    
    UserInput &userInput;
    InSequences inSequences; // when we read a reference we can store it here
    
    uint32_t processedBuffers = 0; // useful to keep track of buffers as they are processed
    uint8_t k; // klen
    uint64_t totKmers = 0, totKmersUnique = 0, totKmersDistinct = 0; // summary statistics
    std::atomic<bool> readingDone{false};
    std::vector<std::thread> threads;
    std::vector<std::future<bool>> futures;
    std::string DBextension;
    
    const uint16_t mapCount = 128; // number of maps to store the kmers, the longer the kmers, the higher number of maps to increase efficiency
    
    const uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount; // this value allows to assign any kmer to a map based on its hashed value
    
    uint64_t* pows = new uint64_t[k]; // storing precomputed values of each power significantly speeds up hashing

    std::vector<Buf<uint8_t>*> buffersVec; // a vector for all buffers
    
    using parallelMap = phmap::parallel_flat_hash_map<uint64_t, TYPE1,
                                              std::hash<uint64_t>,
                                              std::equal_to<uint64_t>,
                                              std::allocator<std::pair<const uint64_t, TYPE1>>,
                                              8,
                                              phmap::NullMutex>;
    
    using parallelMap32 = phmap::parallel_flat_hash_map<uint64_t, TYPE2,
                                              std::hash<uint64_t>,
                                              std::equal_to<uint64_t>,
                                              std::allocator<std::pair<const uint64_t, TYPE2>>,
                                              8,
                                              phmap::NullMutex>;
    
    std::vector<parallelMap*> maps; // all hash maps where TYPE1S are stored
    std::vector<parallelMap32*> maps32;
    
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
    
    std::mutex readMtx, hashMtx;
    
public:
    
    Kmap(UserInput& userInput) : userInput(userInput), k{userInput.kmerLen} {
        
        DBextension = "kc";
        
        for(uint8_t p = 0; p<k; ++p) // precomputes the powers of k
            pows[p] = (uint64_t) pow(4,p);
        
        for(uint16_t m = 0; m<mapCount; ++m)
            maps.push_back(new parallelMap);
        
        for(uint16_t m = 0; m<mapCount; ++m)
            maps32.push_back(new parallelMap32);
        
        if (userInput.kmerDB.size() == 0) { // if we are not reading an existing db
            lg.verbose("Deleting any tmp file");
            for(uint16_t m = 0; m<mapCount; ++m) {// remove tmp buffers and maps if any
                threadPool.queueJob([=]{ return remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str()); });
                threadPool.queueJob([=]{ return remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str()); });
                uint8_t fileNum = 0;
                while (fileExists(userInput.prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum++) +  ".tmp.bin"))
                    threadPool.queueJob([=]{ return remove((userInput.prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin").c_str()); });
                remove((userInput.prefix + "/.index").c_str());
            }
            jobWait(threadPool);
            initHashing(); // start parallel hashing
        }
        
    };
    
    ~Kmap(){ // always need to call the destructor and delete for any object called with new to avoid memory leaks
        for (parallelMap* map : maps)
            delete map;
        for (parallelMap32* map : maps32)
            delete map;
        delete[] pows;
        
    }
    
    uint64_t mapSize(parallelMap& m);
    
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
    
    bool mergeSubMaps(parallelMap* map1, parallelMap* map2, uint8_t subMapIndex, uint16_t m);
    
    bool unionSum(parallelMap* map1, parallelMap* map2, uint16_t m);
    
    bool traverseInReads(std::string* readBatch);
    
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


template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::memoryOk() {
    
    return get_mem_inuse(3) < maxMem;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::memoryOk(int64_t delta) {
    
    return get_mem_inuse(3) + convert_memory(delta, 3) < maxMem;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
uint64_t Kmap<INPUT, TYPE1, TYPE2>::mapSize(parallelMap& m) {
    
   return m.capacity() * (sizeof(typename parallelMap::value_type) + 1) + sizeof(parallelMap);
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::joinThreads() {
    
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

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::initHashing(){
    
    std::packaged_task<bool()> task([this] { return dumpBuffers(); });
    futures.push_back(task.get_future());
    threads.push_back(std::thread(std::move(task)));
    
    int16_t threadN = threadPool.totalThreads() - 1; // substract the writing thread
    
    if (threadN == 0)
        threadN = 1;
    
    for (uint8_t t = 0; t < threadN; t++) {
        std::packaged_task<bool()> task([this] { return hashSequences(); });
        futures.push_back(task.get_future());
        threads.push_back(std::thread(std::move(task)));
    }
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::dumpBuffers() {
    
    bool hashing = true;
    std::vector<Buf<uint8_t>*> buffersVecCpy;
    std::ofstream bufFile[mapCount];
    
    for (uint16_t b = 0; b<mapCount; ++b) // we open all files at once
        bufFile[b] = std::ofstream(userInput.prefix + "/.buf." + std::to_string(b) + ".bin", std::fstream::app | std::ios::out | std::ios::binary);
    
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
        
        for (Buf<uint8_t>* buffers : buffersVecCpy) { // for each array of buffers
            
            for (uint16_t b = 0; b<mapCount; ++b) { // for each buffer file
                
                Buf<uint8_t>* buffer = &buffers[b];
                bufFile[b].write(reinterpret_cast<const char *>(&buffer->pos), sizeof(uint64_t));
                bufFile[b].write(reinterpret_cast<const char *>(buffer->seq), sizeof(uint8_t) * buffer->pos);
                delete[] buffers[b].seq;
                freed += buffers[b].size * sizeof(uint8_t);
                
            }
            delete[] buffers;
        }
    }
    
    for (uint16_t b = 0; b<mapCount; ++b) // we close all files
        bufFile[b].close();
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::buffersToMaps() {
    
    std::vector<std::function<bool()>> jobs;
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of buf files
        fileSizes.push_back(fileSize(userInput.prefix + "/.buf." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        jobs.push_back([this, i] { return processBuffers(i); });
        
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
    
    consolidateTmpMaps();
    dumpHighCopyKmers();

}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::processBuffers(uint16_t m) {
    
    uint64_t pos = 0, hash;
    Buf<uint8_t> *buf;
    
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
        
        parallelMap& map = *maps[m]; // the map associated to this buffer
        parallelMap32& map32 = *maps32[m];
        uint64_t map_size = mapSize(map);
        
        bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
        
        buf = new Buf<uint8_t>(pos);
        
        buf->pos = pos;
        buf->size = pos;
        
        bufFile.read(reinterpret_cast<char *>(buf->seq), sizeof(uint8_t) * buf->pos);
        
        for (uint64_t c = 0; c<pos; c+=8) {
            
            memcpy(&hash, &buf->seq[c], 8);
            
            uint8_t &count = map[hash];
            bool overflow = (count >= 254 ? true : false);
            
            if (!overflow)
                ++count; // increase kmer coverage
            
            if (overflow) {
                
                uint32_t &count32 = map32[hash];
                
                if (count32 == 0) { // first time we add the kmer
                    count32 = count;
                    count = 255; // invalidates int8 kmer
                }

                if (count32 < LARGEST)
                    ++count32; // increase kmer coverage
            }
        }
        
        delete[] buf->seq;
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

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::mergeTmpMaps(uint16_t m) { // a single job merging maps with the same hashes
    
    std::string prefix = userInput.prefix; // loads the first map
    std::string firstFile = prefix + "/.map." + std::to_string(m) + ".0.tmp.bin";
    
    if (!fileExists(prefix + "/.map." + std::to_string(m) + ".1.tmp.bin")) {
        
        std::rename(firstFile.c_str(), (prefix + "/.map." + std::to_string(m) + ".bin").c_str());
        return true;
        
    }
    
    uint8_t fileNum = 0;
    
    while (fileExists(prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) + ".tmp.bin")) { // for additional map loads the map and merges it
        
        std::string nextFile = prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum++) + ".tmp.bin"; // loads the next map
        parallelMap* nextMap = new parallelMap;
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

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::reloadMap32(uint16_t m) {
    
    parallelMap& map = *maps[m]; // the map associated to this buffer
    parallelMap32& map32 = *maps32[m];
    
    for (auto pair : map32) {
        
        uint8_t count = 255;
        auto newPair = std::make_pair(pair.first, count);
        map.insert(newPair);
    }
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::consolidateTmpMaps(){ // concurrent merging of the maps that store the same hashes
    
    lg.verbose("Consolidating temporary maps");
    
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
        fileSizes.push_back(fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".0.tmp.bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        mergeTmpMaps(i);
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::dumpTmpMap(std::string prefix, uint16_t m) {
    
    uint8_t fileNum = 0;
    
    while (fileExists(prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin"))
        ++fileNum;
        
    prefix.append("/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m]->phmap_dump(ar_out);
    
    deleteMap(m);
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::deleteMap(uint16_t m) {
    
    uint64_t map_size = mapSize(*maps[m]);
    delete maps[m];
    freed += map_size;
    
    maps[m] = new parallelMap;
    alloc += mapSize(*maps[m]);
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::dumpHighCopyKmers() {
    
    parallelMap32 map32Total;
    
    for (uint16_t m = 0; m<mapCount; ++m) {
        for (auto pair : *maps32[m])
            map32Total.insert(pair);
        delete maps32[m];
        maps32[m] = new parallelMap32;
    }

    phmap::BinaryOutputArchive ar_out((userInput.prefix + "/.map.hc.bin").c_str());
    map32Total.phmap_dump(ar_out);
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::status() {
    
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - past;
    
    if (elapsed.count() > 0.1) {
        lg.verbose("Read batches: " + std::to_string(readBatches.size()) + ". Hash buffers: " + std::to_string(buffersVec.size()) + ". Memory in use/allocated/total: " + std::to_string(get_mem_inuse(3)) + "/" + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
    
        past = std::chrono::high_resolution_clock::now();
    }
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::loadMap(std::string prefix, uint16_t m) { // loads a specific map
    
    prefix.append("/.map." + std::to_string(m) + ".bin");
    phmap::BinaryInputArchive ar_in(prefix.c_str());
//    allocMemory(fileSize(prefix));
    maps[m]->phmap_load(ar_in);
    alloc += mapSize(*maps[m]);
    
    return true;

}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::loadHighCopyKmers() {
    
    parallelMap32 map32Total;
    phmap::BinaryInputArchive ar_in((userInput.prefix + "/.map.hc.bin").c_str());
    map32Total.phmap_load(ar_in);
    
    for (auto pair : map32Total) {
        uint64_t i = pair.first % mapCount;
        maps32[i]->insert(pair);
    }
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::kunion(){ // concurrent merging of the maps that store the same hashes
    
    parallelMap32 map32Total; // first merge high-copy kmers
    
    for (unsigned int i = 0; i < userInput.kmerDB.size(); ++i) { // for each kmerdb loads the map and merges it
        
        std::string prefix = userInput.kmerDB[i]; // loads the next map
        prefix.append("/.map.hc.bin");
        
        parallelMap32 nextMap;
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
        uint64_t i = pair.first % mapCount;
        maps32[i]->insert(pair);
    }
    
    std::vector<std::function<bool()>> jobs;
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
        fileSizes.push_back(fileSize(userInput.kmerDB[0] + "/.map." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        mergeMaps(i);
    
    dumpHighCopyKmers();
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::mergeMaps(uint16_t m) { // a single job merging maps with the same hashes
    
    std::string prefix = userInput.kmerDB[0]; // loads the first map
    prefix.append("/.map." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    maps[m]->phmap_load(ar_in);

    for (unsigned int i = 1; i < userInput.kmerDB.size(); ++i) { // for each kmerdb loads the map and merges it
        
        std::string prefix = userInput.kmerDB[i]; // loads the next map
        prefix.append("/.map." + std::to_string(m) + ".bin");
        
        parallelMap* nextMap = new parallelMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap->phmap_load(ar_in);
        
        unionSum(nextMap, maps[m], m); // unionSum operation between the existing map and the next map
        delete nextMap;
        
    }
    
    dumpMap(userInput.prefix, m);
    deleteMap(m);
    
    summary(m);
    
    return true;

}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.map." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m]->phmap_dump(ar_out);
    
    deleteMap(m);
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::report() { // generates the output from the program
    
    const static phmap::parallel_flat_hash_map<std::string,int> string_to_case{ // different outputs available
        {"hist",1},
        {"kc",2}

    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    stats();
    
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
            *ostream<<+k<<"\n"<<mapCount<<std::endl;
            ofs.close();
            break;
        }
        default: {}
    }
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::traverseInReads(std::string* readBatch) { // specialized for string objects
    
    while(freeMemory) {status();}
    
    {
        std::lock_guard<std::mutex> lck(readMtx);
        readBatches.push(readBatch);
    }
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
inline uint64_t Kmap<INPUT, TYPE1, TYPE2>::hash(uint8_t *kmer, bool *isFw) { // hashing function for kmers
    
    uint64_t fw = 0, rv = 0; // hashes for both forward and reverse complement sequence
    
    for(uint8_t c = 0; c<k; ++c) { // for each position up to klen
        fw += *kmer * pows[c]; // base * 2^N
        rv += (3-(*kmer++)) * pows[k-c-1]; // we walk the kmer backward to compute the rvcp
    }
    
    if (isFw != NULL)
        *isFw = fw < rv ? true : false; // we preserve the actual orientation for DBG applications
    
    return fw < rv ? fw : rv;
}

template<class INPUT, typename TYPE1, typename TYPE2>
inline std::string Kmap<INPUT, TYPE1, TYPE2>::reverseHash(uint64_t hash) { // hashing function for kmers
    
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

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::consolidate() { // to reduce memory footprint we consolidate the buffers as we go
    status();
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::finalize() { // ensure we count all residual buffers
    if (userInput.kmerDB.size() == 0) {
        readingDone = true;
        joinThreads();
        lg.verbose("Converting buffers to maps");
        buffersToMaps();
    }
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::stats() {
    
    lg.verbose("Computing summary statistics");
    
    std::vector<std::function<bool()>> jobs;
    std::array<uint16_t, 2> mapRange = {0,0};
    
    while (mapRange[1] < mapCount) {
        
        mapRange = computeMapRange(mapRange);
        loadMapRange(mapRange);
        
        for (uint32_t i = mapRange[0]; i < mapRange[1]; ++i)
            jobs.push_back([this, i] { return summary(i); });
        
        threadPool.queueJobs(jobs);
        jobWait(threadPool);
        jobs.clear();
        deleteMapRange(mapRange);
    }
    DBstats();
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::DBstats() {
    
    uint64_t missing = pow(4,k)-totKmersDistinct;
    
    std::cout<<"DB Summary statistics:\n"
             <<"Total kmers: "<<totKmers<<"\n"
             <<"Unique kmers: "<<totKmersUnique<<"\n"
             <<"Distinct kmers: "<<totKmersDistinct<<"\n"
             <<"Missing kmers: "<<missing<<"\n";
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::summary(uint16_t m) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0;
    phmap::parallel_flat_hash_map<uint64_t, uint64_t> hist;
    
    for (auto pair : *maps[m]) {
        
        if (pair.second == 255) // check the large table
            continue;
        
        if (pair.second == 1)
            ++kmersUnique;
        
        ++kmersDistinct;
        ++hist[pair.second];
    }
    
    for (auto pair : *maps32[m]) {
        
        ++kmersDistinct;
        ++hist[pair.second];
    }
 
    std::lock_guard<std::mutex> lck(mtx);
    totKmersUnique += kmersUnique;
    totKmersDistinct += kmersDistinct;
    
    for (auto pair : hist) {
        
        finalHistogram[pair.first] += pair.second;
        totKmers += pair.first * pair.second;
    }
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::printHist(std::unique_ptr<std::ostream>& ostream) { // prints the histogram
    
    std::vector<std::pair<uint64_t, uint64_t>> table(finalHistogram.begin(), finalHistogram.end()); // converts the hashmap to a table
    std::sort(table.begin(), table.end());
    
    for (auto pair : table)
        *ostream<<pair.first<<"\t"<<pair.second<<"\n";

}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::hashSequences() {
    //   Log threadLog;
    std::string *readBatch;
    uint64_t len;
    
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
        }

        if (len<k) {
            delete readBatch;
            freed += len * sizeof(char);
            continue;
        }
        
        Buf<uint8_t> *buffers = new Buf<uint8_t>[mapCount];
        const unsigned char *first = (unsigned char*) readBatch->c_str();
        allocMemory(len * sizeof(uint8_t));
        uint8_t *str = new uint8_t[len];
        uint8_t e = 0;
        uint64_t key, pos = 0, kcount = len-k+1;
        bool isFw = false;
        Buf<uint8_t>* buffer;
        
        for (uint64_t p = 0; p<kcount; ++p) {
            
            for (uint8_t c = e; c<k; ++c) { // generate k bases if e=0 or the next if e=k-1
                
                str[p+c] = ctoi[*(first+p+c)]; // convert the next base
                if (str[p+c] > 3) { // if non-canonical base is found
                    p = p+c; // move position
                    e = 0; // reset base counter
                    break;
                }
                e = k-1;
            }
            
            if (e == 0) // not enough bases for a kmer
                continue;

            key = hash(str+p, &isFw);

            buffer = &buffers[key % mapCount];
            pos = buffer->newPos(8);
            memcpy(&buffer->seq[pos-8], &key, 8);
        }
        delete[] str;
        delete readBatch;
        //    threadLog.add("Processed sequence: " + sequence->header);
        //    std::lock_guard<std::mutex> lck(mtx);
        //    logs.push_back(threadLog);
        std::lock_guard<std::mutex> lck(hashMtx);
        freed += len * sizeof(char) * 2;
        buffersVec.push_back(buffers);
    }
    return true;
}

template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::mergeSubMaps(parallelMap* map1, parallelMap* map2, uint8_t subMapIndex, uint16_t m) {
    
    auto& inner = map1->get_inner(subMapIndex);   // to retrieve the submap at given index
    auto& submap1 = inner.set_;        // can be a set or a map, depending on the type of map1
    auto& inner2 = map2->get_inner(subMapIndex);
    auto& submap2 = inner2.set_;
    parallelMap32& map32 = *maps32[m];
    
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


template<class INPUT, typename TYPE1, typename TYPE2>
bool Kmap<INPUT, TYPE1, TYPE2>::unionSum(parallelMap* map1, parallelMap* map2, uint16_t m) {
    
    std::vector<std::function<bool()>> jobs;
    
    if (map1->subcnt() != map2->subcnt()) {
        fprintf(stderr, "Maps don't have the same numbers of submaps (%zu != %zu). Terminating.\n", map1->subcnt(), map2->subcnt());
        exit(EXIT_FAILURE);
    }
    
    for(std::size_t subMapIndex = 0; subMapIndex < map1->subcnt(); ++subMapIndex)
        jobs.push_back([this, map1, map2, subMapIndex, m] { return mergeSubMaps(map1, map2, subMapIndex, m); });
    
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
    
    return true;
    
}

template<class INPUT, typename TYPE1, typename TYPE2>
std::array<uint16_t, 2> Kmap<INPUT, TYPE1, TYPE2>::computeMapRange(std::array<uint16_t, 2> mapRange) {
    
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

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::loadMapRange(std::array<uint16_t, 2> mapRange) {
    
    std::vector<std::function<bool()>> jobs;
    
    for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
        jobs.push_back([this, m] { return loadMap(userInput.prefix, m); });
    
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::deleteMapRange(std::array<uint16_t, 2> mapRange) {
    
    for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
        deleteMap(m);
}

template<class INPUT, typename TYPE1, typename TYPE2>
void Kmap<INPUT, TYPE1, TYPE2>::cleanup() {
    
    if(!(userInput.kmerDB.size() == 1) && userInput.outFile.find("." + DBextension) == std::string::npos) {
        
        lg.verbose("Deleting tmp files");
        
        std::vector<std::function<bool()>> jobs;
        
        for(uint16_t m = 0; m<mapCount; ++m) // remove tmp files
            jobs.push_back([this, m] { return remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str()); });
        
        threadPool.queueJobs(jobs);
        
        remove((userInput.prefix + "/.map.hc.bin").c_str());
        
        if (userInput.prefix != ".")
            rm_dir(userInput.prefix.c_str());
        
    }
    
    jobWait(threadPool);
    
}

#endif //KMER
