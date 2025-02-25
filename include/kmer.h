#ifndef KMER_H
#define KMER_H

#include "parallel-hashmap/phmap.h"
#include "parallel-hashmap/phmap_dump.h"

#include <fastx.h>
#include <fcntl.h>
#include <MinScan.h>

#define LARGEST 4294967295 // 2^32-1

template<typename T>
inline void freeContainer(T& p_container) // this is a C++ trick to empty a container and release associated memory
{
	T empty;
	std::swap(p_container, empty); // swapping a container with an empty (NULL) container should release associated memory
}

class Key {
	uint64 offset;
public:
	
	Key(){}
	Key(uint64 offset) : offset(offset){}
	
	void assignOffset(uint64 offset) {
		this->offset = offset;
	}
	uint64 getOffset() const {
		return offset;
	}
};

struct HcKmer {
	uint64_t key;
	uint32_t value;
	uint8_t map;
};

struct KeyHasher {

	uint64 *data;
	
	KeyHasher() {}
	
	KeyHasher(uint64 *data) : data(data) {}
	
	std::size_t operator()(const Key& key) const {

		uint64 hash;
		Get_Hash(&hash, data, key.getOffset());
		return hash;
	}
};

struct KeyEqualTo {
	
	uint32_t k;
	uint64 *data;
	
	KeyEqualTo() {}
	
	KeyEqualTo(uint64 *data, uint32_t k) : k(k), data(data) {}
	
	bool operator()(const Key& key1, const Key& key2) const {
		
		uint64 hash1, hash2;
		int dir1 = Get_Hash(&hash1, data, key1.getOffset());
		int dir2 = Get_Hash(&hash2, data, key2.getOffset());
		
		if (k <= 32)
			return hash1 == hash2;
		else if (hash1 != hash2)
			return false;
		
		uint64 *kmer1 = New_Supermer_Buffer(), *kmer2 = New_Supermer_Buffer();
		
		Get_Canonical_Kmer(kmer1,dir1,hash1,data,key1.getOffset());
		Get_Canonical_Kmer(kmer2,dir2,hash2,data,key2.getOffset());
		
		int w = 0;
		int x = 62;
		for (uint32_t i = 33; i < k; i++) {
			
		if (((kmer1[w]>>x)&0x3llu) != ((kmer2[w]>>x)&0x3llu))
				return false;
			
			if (x == 0)
				{ w += 1;
				  x = 64;
			  }
			x -= 2;
		}
		
		delete kmer1;
		delete kmer2;
		
		return true;
	}
};

struct SeqBuf {
	uint64 *data;
	uint64_t len;
};

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2> // DERIVED implements the CRTP technique, INPUT is a specialized userInput type depending on the tool, KEY is the key type for the hashtable, TYPE1 is the low frequency type of elements we wish to store in the maps, e.g. uint8_t kmer counts, TYPE2 is the high frequency type of elements we wish to store in the maps, e.g. uint32_t kmer counts
class Kmap {

protected: // they are protected, so that they can be further specialized by inheritance
	
	UserInput &userInput;
	InSequences inSequences; // when we read a reference we can store it here
	
	uint8_t sLen; // smer length
	uint32_t k; // kmer length
	uint64_t tot = 0, totUnique = 0, totDistinct = 0; // summary statistics
	std::atomic<bool> readingDone{false};
	std::string DBextension;
	
	const static uint16_t mapCount = 128; // number of maps to store the kmers, the longer the kmers, the higher number of maps to increase efficiency
	
	const uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount; // this value allows to assign any kmer to a map based on its hashed value
	
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
	
	std::vector<ParallelMap*> maps; // all hash maps where TYPE1 are stored
	std::vector<ParallelMap32*> maps32; // all hash maps where TYPE2 are stored
	
	std::vector<ParallelMap*> tmpMaps[mapCount];
	std::vector<ParallelMap32*> tmpMaps32[mapCount];
	
	std::vector<bool> mapsInUse = std::vector<bool>(mapCount, false); // useful with multithreading to ensure non-concomitant write access to maps
	
	phmap::flat_hash_map<uint64_t, uint64_t> finalHistogram; // the final kmer histogram
	
	const unsigned char ctoi[256] = { // this converts ACGT>0123 and any other character to 4 in time O(1)
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
	
	uint8_t mapReady = 0, toDelete = 0, deleted = 0;
	std::array<uint16_t,mapCount> hashBufferDone{}, mapDoneCounts{};
	
	std::mutex readMtx, hashMtx, summaryMtx;
	std::condition_variable readMutexCondition, hashMutexCondition, summaryMtxCondition;
	uint16_t bufferDone = 0;
	
	int bufferFiles[mapCount];
	SeqBuf seqBuf[mapCount];
	
public:
	
	Kmap(UserInput& userInput) : userInput(userInput), sLen{userInput.sLen}, k{userInput.kLen} {
		
		DBextension = "kc";
		
		maps.reserve(mapCount);
		maps32.reserve(mapCount);
		
		if (userInput.kmerDB.size() == 0) { // if we are not reading an existing db
			lg.verbose("Deleting any tmp file");
			for(uint16_t m = 0; m<mapCount; ++m) {// remove tmp buffers and maps if any
				remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str());
				remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str());
				uint8_t fileNum = 0;
				while (fileExists(userInput.prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin")) {
					remove((userInput.prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin").c_str());
					++fileNum;
				}
				remove((userInput.prefix + "/.index").c_str());
				remove((userInput.prefix + "/.seq.bin").c_str());

			}
			for(uint16_t m = 0; m<mapCount; ++m)
				bufferFiles[m] = open((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str(), O_WRONLY | O_CREAT, S_IRWXU | S_IRWXG | S_IRWXO);
			
			initBuffering(); // start parallel buffering
		}
	};
	
	uint64_t mapSize(ParallelMap& m);

	
	bool memoryOk();
	
	bool memoryOk(int64_t delta);
	
	void initBuffering();
	
	void buffersToMaps();
	
	void initHashing();
	
	bool hashBuffer(uint16_t thisThread);
	
	bool consolidateTmpMap(uint16_t m);
	
	bool dumpTmpMap(std::string prefix, uint8_t m, ParallelMap *map);
	
	bool deleteMap(uint16_t m);
	
	void dumpHighCopyKmers();
	
	bool mergeTmpMaps(uint16_t m);
	
	void status();
	
	void status2(uint8_t buffers);
	
	void kunion();
	
	bool mergeSubMaps(ParallelMap* map1, ParallelMap* map2, uint8_t subMapIndex, uint16_t m);
	
	bool unionSum(ParallelMap* map1, ParallelMap* map2, uint16_t m);
	
	bool traverseInReads(std::string* readBatch);
	
	bool traverseInReads(Sequences* readBatch);
	
	inline uint64_t hash(Buf2bit<> *kmerPtr, uint64_t p, bool* isFw = NULL);
	
	inline uint64_t hash(uint64 *kmerPtr, uint32_t p, bool* isFw = NULL);
	
	inline std::string reverseHash(uint64_t hash);
	
	bool generateBuffers();
	
	void consolidate();
	
	void finalize();
	
	void computeStats();
	
	void DBstats();
	
	bool summary(uint16_t m);
	
	void printHist(std::unique_ptr<std::ostream>& ostream);
	
	void report();
	
	void loadBuffer(uint8_t m);
	
	bool dumpMap(std::string prefix, uint16_t m);
	
	bool loadMap(std::string prefix, uint16_t m);
	
	std::array<uint16_t, 2> computeMapRange(std::array<uint16_t, 2> mapRange);
	
	void loadMapRange(std::array<uint16_t, 2> mapRange);
	
	void deleteMapRange(std::array<uint16_t, 2> mapRange);
	
	void cleanup();
	
	bool mergeMaps(uint16_t m);
	
	void mergeMaps(ParallelMap &map1, ParallelMap &map2, ParallelMap32 &map32);

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
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::generateBuffers() {
	//   Log threadLog;
	std::string *readBatch;
	uint64_t len = 0;

	Distribution_Bundle* bundle = Begin_Distribution(bufferFiles);
	
	while (true) {
			
		{
			std::unique_lock<std::mutex> lck(readMtx);
			
			if (readingDone && readBatches.size() == 0) {
				End_Distribution(bundle);
				++bufferDone;
				readMutexCondition.notify_all();
				return true;
			}

			if (readBatches.size() == 0)
				continue;
			
			readBatch = readBatches.front();
			readBatches.pop();
			len = readBatch->size();
		}
		
		Distribute_Sequence(const_cast<char*>(readBatch->data()), len, bundle);
		
		delete readBatch;
		freed += len * sizeof(char);
	}
	return true;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::initBuffering(){
	
	Init_Genes_Package(k, sLen);
	
	std::vector<std::function<bool()>> jobs;
	for (uint8_t t = 0; t < threadPool.totalThreads(); ++t)
		jobs.push_back([this] { return static_cast<DERIVED*>(this)->generateBuffers(); });
	
	threadPool.queueJobs(jobs);
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::initHashing(){
	
	std::vector<std::function<bool()>> jobs;
	for (uint8_t t = 0; t < threadPool.totalThreads(); ++t)
		jobs.push_back([this, t] { return static_cast<DERIVED*>(this)->hashBuffer(t); });
	threadPool.queueJobs(jobs);
	
	for (uint16_t m = 0; m<mapCount; ++m) { // the master thread reads the buffers in
		
		loadBuffer(m);
		tmpMaps[m].resize(threadPool.totalThreads());
		tmpMaps32[m].resize(threadPool.totalThreads());
		{
			std::lock_guard<std::mutex> lck(hashMtx);
			++mapReady;
		}
		hashMutexCondition.notify_all();
		{
			std::unique_lock<std::mutex> lck(summaryMtx);
			summaryMtxCondition.wait(lck, [&] {
				if (mapDoneCounts[m] == threadPool.totalThreads())
					return true;
				return false;
			});
		}
		threadPool.queueJob([=]{ return consolidateTmpMap(m); });
	}
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::hashBuffer(uint16_t t) {
	
	int count;
	uint64 *data;
	uint64 offset, hash;
	uint32_t totalThreads = threadPool.totalThreads();
	uint8_t m = 0;
	
	while (m < mapCount) {
		{
			std::unique_lock<std::mutex> lck(hashMtx);
			hashMutexCondition.wait(lck, [this,m] {
				return mapReady > m;
			});
		}
		data = seqBuf[m].data;

		Scan_Bundle *bundle = Begin_Supermer_Scan(data, seqBuf[m].len);
		uint64 *super = New_Supermer_Buffer();
		
		tmpMaps[m][t] = new ParallelMap(0, KeyHasher(data), KeyEqualTo(data, k));
		tmpMaps[m][t]->reserve(seqBuf[m].len);
		tmpMaps32[m][t] = new ParallelMap32(0, KeyHasher(data), KeyEqualTo(data, k));
		
		ParallelMap &map = *tmpMaps[m][t];
		ParallelMap32 &map32 = *tmpMaps32[m][t];

		while ((count = Get_Kmer_Count(bundle))) {
			
			offset = Current_Offset(data, bundle);
			
			for (int c = 0; c<count; ++c) {
				
				Get_Hash(&hash, data, offset);
				if (hash % totalThreads != t)
					continue;
					
				Key key(offset);
				TYPE1 &count = map[key];
				bool overflow = (count >= 254 ? true : false);
				
				if (!overflow)
					++count; // increase kmer coverage
				else {
					
					TYPE2 &count32 = map32[key];
					
					if (count32 == 0) { // first time we add the kmer
						count32 = count;
						count = 255; // invalidates int8 kmer
					}
					if (count32 < LARGEST)
						++count32; // increase kmer coverage
				}
				offset += 2;
			}
			Skip_Kmers(count, bundle);
		}
		free(super);
		End_Supermer_Scan(bundle);
		{
			std::lock_guard<std::mutex> lck(summaryMtx);
			++mapDoneCounts[m];
		}
		summaryMtxCondition.notify_one();
		++m;
	}
	return true;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::loadBuffer(uint8_t m){
	uint64_t pos;
	
	std::string fl = userInput.prefix + "/.buf." + std::to_string(m) + ".bin";
	if (!fileExists(fl)) {
		fprintf(stderr, "Buffer file %s does not exist. Terminating.\n", fl.c_str());
		exit(EXIT_FAILURE);
	}
	std::ifstream bufFile(fl, std::ios::binary | std::ios::ate);
	std::streamsize size = bufFile.tellg();
	bufFile.seekg(0, std::ios::beg);
	pos = size/8;
	
	uint64 *data = new uint64[pos];
	if (!bufFile.read(reinterpret_cast<char *>(data),  sizeof(uint8_t) * size)) {
			printf("Cannot read %i buffer. Terminating.\n", m);
		exit(EXIT_FAILURE);
	}
	bufFile.close();
	
	seqBuf[m].data = data;
	seqBuf[m].len = pos;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::buffersToMaps() {
	
	initHashing();
	jobWait(threadPool);

	//consolidateTmpMaps();
	dumpHighCopyKmers();
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
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::consolidateTmpMap(uint16_t m){ // concurrent merging of the maps that store the same hashes
	
	maps[m] = new ParallelMap(0, KeyHasher(seqBuf[m].data), KeyEqualTo(seqBuf[m].data, k));
	maps[m]->reserve(seqBuf[m].len);
	maps32[m] = new ParallelMap32(0, KeyHasher(seqBuf[m].data), KeyEqualTo(seqBuf[m].data, k));
	
	for (uint32_t t = 0; t < tmpMaps[m].size(); ++t) {
		maps[m]->insert(tmpMaps[m][t]->begin(), tmpMaps[m][t]->end());
		delete tmpMaps[m][t];
		maps32[m]->insert(tmpMaps32[m][t]->begin(), tmpMaps32[m][t]->end());
		delete tmpMaps32[m][t];
	}
	summary(m);
	delete seqBuf[m].data;
	dumpTmpMap(userInput.prefix, m, maps[m]);
	return true;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::dumpTmpMap(std::string prefix, uint8_t m, ParallelMap *map) {
	
	uint8_t fileNum = 0;
	
	while (fileExists(prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin"))
		++fileNum;
		
	prefix.append("/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin");
	
	phmap::BinaryOutputArchive ar_out(prefix.c_str());
	
	if (map == NULL) {
		map = new ParallelMap;
		map->phmap_dump(ar_out);
		delete map;
		map = NULL;
	}else{
		map->phmap_dump(ar_out);
		delete map;
	}
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
	
	std::ofstream bufFile = std::ofstream(userInput.prefix + "/.hc.bin", std::ios::out | std::ios::binary);
	
	for (uint16_t m = 0; m<mapCount; ++m) { // write number of hc kmers for each map
		uint64_t count = maps32[m]->size();
		bufFile.write(reinterpret_cast<const char *>(&count), sizeof(uint64_t));
	}
		
	for (uint16_t m = 0; m<mapCount; ++m) {

		for (auto pair : *maps32[m]) {
			HcKmer hcKmer;
			hcKmer.key = pair.first.getOffset();
			hcKmer.value = pair.second;
			hcKmer.map = m;
			bufFile.write(reinterpret_cast<const char *>(&hcKmer.key), sizeof(uint64_t));
			bufFile.write(reinterpret_cast<const char *>(&hcKmer.value), sizeof(uint32_t));
			bufFile.write(reinterpret_cast<const char *>(&hcKmer.map), sizeof(uint8_t));
		}
		delete maps32[m];
	}
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::status() {
	
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - past;
	
	if (elapsed.count() > 0.1) {
		lg.verbose("Read batches: " + std::to_string(readBatches.size() + sequenceBatches.size()) + ". Memory in use/allocated/total: " + std::to_string(get_mem_inuse(3)) + "/" + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
	
		past = std::chrono::high_resolution_clock::now();
	}
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::status2(uint8_t buffers) {
	
//    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - past;
	
//    if (elapsed.count() > 0.1) {
		lg.verbose("Buffers: " + std::to_string(buffers) + ". Memory in use/allocated/total: " + std::to_string(get_mem_inuse(3)) + "/" + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
	
//        past = std::chrono::high_resolution_clock::now();
//    }
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::loadMap(std::string prefix, uint16_t m) { // loads a specific map
	
	loadBuffer(m);
	prefix.append("/.map." + std::to_string(m) + ".bin");
	phmap::BinaryInputArchive ar_in(prefix.c_str());
	maps[m] = new ParallelMap(0, KeyHasher(seqBuf[m].data), KeyEqualTo(seqBuf[m].data, k));
	maps[m]->phmap_load(ar_in);
	alloc += mapSize(*maps[m]);
	
	std::ifstream bufFile = std::ifstream(userInput.prefix + "/.hc.bin", std::ios::in | std::ios::binary);
	
	uint64_t offset = 0, counts = 0;
	for(uint16_t i = 0; i<mapCount; ++i) { // find position for this map
		uint64_t count;
		bufFile.read(reinterpret_cast<char *>(&count), sizeof(uint64_t));
		if (i<m)
			offset += count;
		if(i == m)
			counts = count;
	}

	bufFile.seekg(mapCount*sizeof(uint64_t) + offset*(sizeof(uint64_t)+sizeof(uint32_t)+sizeof(uint8_t)));
	
	maps32[m] = new ParallelMap32(0, KeyHasher(seqBuf[m].data), KeyEqualTo(seqBuf[m].data, k));
	
	for (uint64_t i = 0; i<counts; ++i) { // load respective hc kmers
		HcKmer hcKmer;
		bufFile.read(reinterpret_cast<char *>(&hcKmer.key), sizeof(uint64_t));
		bufFile.read(reinterpret_cast<char *>(&hcKmer.value), sizeof(uint32_t));
		bufFile.read(reinterpret_cast<char *>(&hcKmer.map), sizeof(uint8_t));
//        std::cout<<+hcKmer.key<<" "<<+hcKmer.value<<" "<<+hcKmer.map<<std::endl;
		maps32[hcKmer.map]->emplace(std::make_pair(Key(hcKmer.key),hcKmer.value));
	}
	return true;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::kunion(){ // concurrent merging of the maps that store the same hashes
	
	ParallelMap32 map32Total; // first merge high-copy kmers
	
	for (unsigned int i = 0; i < userInput.kmerDB.size(); ++i) { // for each kmerdb loads the map and merges it
		
		std::string prefix = userInput.kmerDB[i]; // loads the next map
		prefix.append("/.hc.bin");
		
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
	
//    for (auto pair : map32Total) {
//        uint64_t pos = pair.first.getOffset();
//        uint64_t i = hash(seqBuf, pos) % mapCount;
//        maps32[i]->emplace(pair);
//    }
	
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
	
//    std::string prefix = userInput.kmerDB[0]; // loads the first map
//    prefix.append("/.map." + std::to_string(m) + ".bin");
//
//    phmap::BinaryInputArchive ar_in(prefix.c_str());
//    maps[m]->phmap_load(ar_in);
//
//    uint64_t pos;
//    std::ifstream bufFile = std::ifstream(userInput.kmerDB[0]+ "/.seq.bin", std::ios::in | std::ios::binary);
//    bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
//    seqBuf[m].seq = new Buf2bit<>(pos);
//    seqBuf[m].seq->size = pos;
//    bufFile.read(reinterpret_cast<char *>(seqBuf[m].seq->seq), sizeof(uint8_t) * seqBuf[m].seq->size);
//
//    for (unsigned int i = 1; i < userInput.kmerDB.size(); ++i) { // for each kmerdb loads the map and merges it
//
//        bufFile = std::ifstream(userInput.kmerDB[0] + "/.seq.bin", std::ios::in | std::ios::binary);
//        bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
//        seqBuf2[m].seq = new Buf2bit<>(pos);
//        seqBuf2[m].seq->size = pos;
//        bufFile.read(reinterpret_cast<char *>(seqBuf2[m].seq->seq), sizeof(uint8_t) * seqBuf2[m].seq->size);
//
//        std::string prefix = userInput.kmerDB[i]; // loads the next map
//        prefix.append("/.map." + std::to_string(m) + ".bin");
//
//        ParallelMap* nextMap = new ParallelMap;
//        phmap::BinaryInputArchive ar_in(prefix.c_str());
//        nextMap->phmap_load(ar_in);
//
//        unionSum(nextMap, maps[m], m); // unionSum operation between the existing map and the next map
//        delete nextMap;
//    }
//
//    dumpMap(userInput.prefix, m);
//    deleteMap(m);
//
//    static_cast<DERIVED*>(this)->summary(m);
	
	return true;

}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::mergeMaps(ParallelMap &map1, ParallelMap &map2, ParallelMap32 &map32) {
	
	for (auto pair : map1) { // for each element in map1, find it in map2 and increase its value
		
		bool overflow = false;
		
		if (pair.second == 255) // already added to int32 map
			continue;
		
		auto got = map32.find(pair.first); // check if this is already a high-copy kmer
		if (got != map32.end()) {
			overflow = true;
		}else{
			TYPE1 &count = map2[pair.first]; // insert or find this kmer in the hash table
				
			if (255 - count <= pair.second)
				overflow = true;
			
			if (!overflow)
				count += pair.second; // increase kmer coverage
		}
		if (overflow) {
			TYPE2 &count32 = map32[pair.first];
			
			if (count32 == 0) { // first time we add the kmer
				TYPE1& count = map2[pair.first];
				count32 = count;
				count = 255; // invalidates int8 kmer
			}
			
			if (LARGEST - count32 >= pair.second)
				count32 += pair.second; // increase kmer coverage
			else
				count32 = LARGEST;
		}
	}
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
		{"stdout",0},
		{"hist",1},
		{"kc",2}

	};
	
	std::string ext = "stdout";
	
	if (userInput.outFile != "")
		ext = getFileExt("." + userInput.outFile);
	
	static_cast<DERIVED*>(this)->DBstats();
	
	lg.verbose("Writing ouput: " + ext);
	
	std::unique_ptr<std::ostream> ostream; // smart pointer to handle any kind of output stream
	
	switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
			
		case 0: { // default
			break;
		}
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
		default: { // unrecognized ext
			std::cout<<"Unrecognized format (."<<ext<<")."<<std::endl;
			exit(EXIT_FAILURE);
		}
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
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::consolidate() { // to reduce memory footprint we consolidate the buffers as we go
	status();
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::finalize() { // ensure we count all residual buffers
	if (userInput.kmerDB.size() == 0) {
		readingDone = true;
		{
			std::unique_lock<std::mutex> lck(readMtx);
			readMutexCondition.wait(lck, [&] {
				status();
				return bufferDone == threadPool.totalThreads();
			});
		}
		lg.verbose("Converting buffers to maps");
		static_cast<DERIVED*>(this)->buffersToMaps();
	}
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::computeStats() {
	
	lg.verbose("Computing summary statistics");
	
	Init_Genes_Package(k, sLen);
	
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
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::DBstats() {
	
	uint64_t missing = pow(4,k)-totDistinct;
	
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
				submap2.emplace(pair);
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
		
		max += fileSize(userInput.prefix + "/.buf." + std::to_string(m) + ".bin");
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
	
	for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m) {
		delete seqBuf[m].data;
		deleteMap(m);
	}
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::cleanup() {
	
	if(!(userInput.kmerDB.size() == 1) && userInput.outFile.find("." + DBextension) == std::string::npos) {
		
		lg.verbose("Deleting tmp files");
		for(uint16_t m = 0; m<mapCount; ++m) { // remove tmp files
			remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str());
			remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str());
		}
		
		remove((userInput.prefix + "/.hc.bin").c_str());
		remove((userInput.prefix + "/.seq.bin").c_str());
		
		if (userInput.prefix != ".")
			rm_dir(userInput.prefix.c_str());
		
	}
	
	jobWait(threadPool);
	
}

#endif //KMER
