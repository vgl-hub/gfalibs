#ifndef KMER_H
#define KMER_H

#include <random>
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
	int shift;
	
	KeyHasher() {}
	
	KeyHasher(uint64 *data, uint32_t k) : data(data), shift((32-k)*2) {
		if (shift < 0)
			shift = 0;
	}
	
	std::size_t operator()(const Key& key) const {

		uint64 hash;
		Get_Hash(&hash, data, key.getOffset());
		return hash >> shift;
	}
};

struct KeyEqualTo {
	
	uint32_t k;
	uint64 *data = nullptr;
	uint64 *kmer1 = nullptr, *kmer2 = nullptr;
	
	KeyEqualTo() {}
	
	KeyEqualTo(uint64 *data, uint32_t k) : k(k), data(data) {
		kmer1 = New_Supermer_Buffer();
		kmer2 = New_Supermer_Buffer();
	}
	
	~KeyEqualTo(){ // need to follow rule of three
//		if (kmer1 != nullptr)
//			free(kmer1);
//		if (kmer2 != nullptr)
//			free(kmer2);
	}
	
	bool operator()(const Key& key1, const Key& key2) const {
		
		uint64 hash1, hash2;
		int dir1 = Get_Hash(&hash1, data, key1.getOffset());
		int dir2 = Get_Hash(&hash2, data, key2.getOffset());
		if (k <= 32)
			return hash1 == hash2;
		else if (hash1 != hash2)
			return false;
		
		Get_Canonical_Kmer(kmer1,dir1,hash1,data,key1.getOffset());
		Get_Canonical_Kmer(kmer2,dir2,hash2,data,key2.getOffset());
		int w = 0;
		int x = 62;
		for (uint32_t i = 0; i < k; ++i) {

			if (((kmer1[w]>>x)&0x3llu) != ((kmer2[w]>>x)&0x3llu))
				return false;

			if (x == 0) {
				++w;
				x = 64;
			}
			x -= 2;
		}
		return true;
	}
};

struct SeqBuf {
	uint64 *data = NULL;
	uint64_t len;
};

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2> // DERIVED implements the CRTP technique, INPUT is a specialized userInput type depending on the tool, KEY is the key type for the hashtable, TYPE1 is the low frequency type of elements we wish to store in the maps, e.g. uint8_t kmer counts, TYPE2 is the high frequency type of elements we wish to store in the maps, e.g. uint32_t kmer counts
class Kmap {

protected: // they are protected, so that they can be further specialized by inheritance
	
	UserInputKcount &userInput;
	InSequences inSequences; // when we read a reference we can store it here
	
	uint8_t sLen; // smer length
	uint32_t k; // kmer length
	uint64_t tot = 0, totUnique = 0, totDistinct = 0; // summary statistics
	std::string DBextension;
	
	const static uint16_t mapCount = 128; // number of maps to store the kmers, the longer the kmers, the higher number of maps to increase efficiency
	
	const uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount; // this value allows to assign any kmer to a map based on its hashed value
	
	using ParallelMap = phmap::flat_hash_map<KEY, TYPE1,
											  KeyHasher,
											  KeyEqualTo,
											  std::allocator<std::pair<const KEY, TYPE1>>>;

	using ParallelMap32 = phmap::flat_hash_map<KEY, TYPE2,
											  KeyHasher,
											  KeyEqualTo,
											  std::allocator<std::pair<const KEY, TYPE2>>>;
	
	std::vector<ParallelMap*> maps; // all hash maps where TYPE1 are stored
	std::vector<ParallelMap32*> maps32; // all hash maps where TYPE2 are stored
	
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
	std::string* readBatches;
	
	uint8_t toDelete = 0, deleted = 0;
	std::array<int16_t,mapCount> mapDoneCounts{};
	std::array<uint64_t,mapCount> mapSizeCounts{};
	
	std::mutex readMtx, hashMtx, summaryMtx;
	std::condition_variable processBufferMutexCondition, hashMutexCondition, summaryMutexCondition;
	
	int bufferFiles[mapCount];
	SeqBuf seqBuf[mapCount];
	
	std::vector<std::thread> writeThreads;
	std::queue<std::pair<uint8_t,uint8_t>> buffersQueue;
	std::queue<uint8_t> finalMapsQueue;
	uint8_t totalMapsDone = 0, bufferProcessed = 0, lastBuffer = 0;
	
public:
	
	Kmap(UserInputKcount& userInput) : userInput(userInput), sLen{userInput.sLen}, k{userInput.kLen} {
		
		DBextension = "kc";
		
		maps.resize(mapCount);
		maps32.resize(mapCount);
		
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
				remove((userInput.prefix + "/.hc.bin").c_str());

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
	
	bool hashBuffer();
	
	bool consolidateTmpMaps();
	
	bool dumpTmpMap(std::string prefix, uint8_t m, ParallelMap *map, uint8_t fileNum);
	
	bool deleteMap(uint16_t m);
	
	void dumpHighCopyKmers();
	
	bool mergeTmpMaps(uint16_t m);
	
	void status();
	
	void status2(uint8_t buffers);
	
	void kunion();
	
	bool mergeSubMaps(ParallelMap* map1, ParallelMap* map2, uint8_t subMapIndex, uint16_t m);
	
	bool unionSum(ParallelMap* map1, ParallelMap* map2, uint16_t m);
	
	bool traverseInReads(std::string* readBatches);
	
	inline uint64_t hash(Buf2bit<> *kmerPtr, uint64_t p, bool* isFw = NULL);
	
	inline uint64_t hash(uint64 *kmerPtr, uint32_t p, bool* isFw = NULL);
	
	inline std::string reverseHash(uint64_t hash);
	
	void readFastqStream(std::shared_ptr<std::istream> input, std::string &buffer);
	
	bool generateBuffers(std::shared_ptr<std::istream> stream);
	
	void getData(std::string *readBatches, std::string &readBatch, uint16_t t);
	
	uint64_t findNextSequence(const std::string *readBatches, const uint64_t quota, const uint16_t t);
	
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
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::initBuffering(){
	
	Init_Genes_Package(k, sLen);
	
	if (userInput.inFiles.size() > 0) {
		
		lg.verbose("Loading input sequences");
		unsigned int numFiles = userInput.inFiles.size();
		
		for (unsigned int i = 0; i < numFiles; ++i) { // for each input read file
			
			//stream objects
			StreamObj streamObj;
			std::shared_ptr<std::istream> stream = streamObj.openStream(userInput, 'r', i);
			
			std::vector<std::function<bool()>> jobs;
			for (uint8_t t = 0; t < threadPool.totalThreads(); ++t)
				jobs.push_back([this, stream] { return static_cast<DERIVED*>(this)->generateBuffers(stream); });
			
			threadPool.queueJobs(jobs);
			jobWait(threadPool);
		}
		lg.verbose("Reads loaded.");
		finalize();
	}else{
		fprintf(stderr, "Reads not provided. Exiting.\n");
		exit(EXIT_FAILURE);
	}
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::consolidate() { // to reduce memory footprint we consolidate the buffers as we go
	status();
}

#define INITIAL_READ_SIZE 1048576  // 1MB
template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::readFastqStream(std::shared_ptr<std::istream> input, std::string &buffer) {

	// Read the first 1MB
	buffer.resize(INITIAL_READ_SIZE);
	input->read(&buffer[0], INITIAL_READ_SIZE);
	std::streamsize bytesRead = input->gcount();
	buffer.resize(bytesRead);  // Resize to actual bytes read
	std::string line;

	while (true) {
	
		std::streampos pos = input->tellg();
		std::getline(*input, line);
		if (!*input)
			break;
		
		// Check if the new line starts a FASTQ record (valid header)
		if (!line.empty() && line[0] == '@') { // it could be a new fastq record
		
			int c = input->peek();  // peek character
			
			if (c == '@') { // this was indeed the end of a quality line
				buffer.insert(buffer.end(), line.begin(), line.end());
				buffer.push_back('\n');
			}else{ // put the line back in the buffer
				input->seekg(pos);
			}
			break;  // Stop reading at the start of the next record
		}
		buffer.insert(buffer.end(), line.begin(), line.end());
		buffer.push_back('\n');
	}
}

#define BUFFER_RESERVE_SIZE 2097152 // 2MB preallocation for efficiency
template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::generateBuffers(std::shared_ptr<std::istream> stream) {
	//   Log threadLog;
	
	std::string readBatch, kmerBatch;
	readBatch.reserve(BUFFER_RESERVE_SIZE); // Preallocate memory for efficiency
	Distribution_Bundle* bundle = Begin_Distribution(bufferFiles);
	
	while (true) {
			
		{
			std::lock_guard<std::mutex> lck(readMtx);
			if (!*stream) {
				End_Distribution(bundle);
				return true;
			}
			readFastqStream(stream, readBatch);
		}
		std::stringstream ss(readBatch);
		getKmers(ss, kmerBatch, readBatch.size());
		Distribute_Sequence(const_cast<char*>(readBatch.data()), readBatch.size(), bundle);
	}
	return true;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::finalize() { // ensure we count all residual buffers
	if (userInput.kmerDB.size() == 0) {
		lg.verbose("Converting buffers to maps");
		static_cast<DERIVED*>(this)->buffersToMaps();
	}
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
void Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::initHashing(){
	
	std::vector<std::function<bool()>> jobs;
	for (uint8_t t = 0; t < threadPool.totalThreads(); ++t)
		jobs.push_back([this] { return static_cast<DERIVED*>(this)->hashBuffer(); });
	threadPool.queueJobs(jobs);
	
	for (uint8_t t = 0; t < userInput.writeThreads; ++t)
		writeThreads.push_back(std::thread(&Kmap::consolidateTmpMaps, this));
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::hashBuffer() {
	
	int count;
	const int hThreads = userInput.hashThreads;
	uint64 *data;
	uint64 offset, hash, mapSizeCount;
	uint8_t m, mapN;
	
	std::vector<uint64_t> fileSizes;
	for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
		fileSizes.push_back(fileSize(userInput.prefix + "/.buf." + std::to_string(m) + ".bin"));
	std::vector<uint32_t> sorted = sortedIndex(fileSizes, true); // sort by largest
	
	while (true) {
		
		{
			std::unique_lock<std::mutex> lck(hashMtx);
			hashMutexCondition.wait(lck, [this,sorted,hThreads] {
				
				if (!buffersQueue.size() && lastBuffer < mapCount) {
					uint8_t idx = lastBuffer++;
					loadBuffer(idx);
					tmpMaps32[idx].resize(hThreads);
					
					for (int32_t i = 0; i < hThreads; ++i)
						buffersQueue.push(std::make_pair(idx, i));
				}
				return (buffersQueue.size() || (bufferProcessed == mapCount));
			});
			if (bufferProcessed == mapCount)
				break;
			std::pair<uint8_t, uint8_t> buf = buffersQueue.front();
			m = buf.first;
			mapN = buf.second;
			//std::cout << +m << " " << +mapN << std::endl;
			buffersQueue.pop();
		}
		hashMutexCondition.notify_all();
		
		data = seqBuf[m].data;
		Scan_Bundle *bundle = Begin_Supermer_Scan(data, seqBuf[m].len);
		
		ParallelMap *ptr = new ParallelMap(0, KeyHasher(data, k), KeyEqualTo(data, k));
		//tmpMaps[m][t]->reserve(seqBuf[m].len);
		tmpMaps32[m][mapN] = new ParallelMap32(0, KeyHasher(data, k), KeyEqualTo(data, k));
		
		ParallelMap &map = *ptr;
		ParallelMap32 &map32 = *tmpMaps32[m][mapN];
		int shift = (32-k)*2;
		if (shift < 0)
			shift = 0;

		while ((count = Get_Kmer_Count(bundle))) {
			
			offset = Current_Offset(data, bundle);
			
			for (int c = 0; c<count; ++c) {
				
				Get_Hash(&hash, data, offset);
				
				if ((hash >> shift) % hThreads == mapN) {
					
					Key key(offset);
					TYPE1 &count = map[key];
					if (count < 254)
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
				}
				offset += 2;
			}
			Skip_Kmers(count, bundle);
		}
		End_Supermer_Scan(bundle);
		mapSizeCount = map.size();
		dumpTmpMap(userInput.prefix, m, &map, mapN);
		{
			std::lock_guard<std::mutex> lck(summaryMtx);
			++mapDoneCounts[m];
			mapSizeCounts[m] += mapSizeCount;
			
			if (mapDoneCounts[m] == hThreads) {
				++bufferProcessed;
				finalMapsQueue.push(m);
				delete seqBuf[m].data;
				summaryMutexCondition.notify_one();
				hashMutexCondition.notify_all();
			}
		}
	}
	summaryMutexCondition.notify_all();
	return true;
}

template<class DERIVED, class INPUT, typename KEY, typename TYPE1, typename TYPE2>
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::consolidateTmpMaps(){ // concurrent merging of the maps that store the same hashes
	
	while (totalMapsDone < mapCount) {
		
		uint8_t m;
		{
			std::unique_lock<std::mutex> lck(summaryMtx);
			summaryMutexCondition.wait(lck, [this] {
				return finalMapsQueue.size() || (totalMapsDone == mapCount);
			});
			if (totalMapsDone == mapCount)
				break;
			m = finalMapsQueue.front();
			finalMapsQueue.pop();
			++totalMapsDone;
		}
		loadBuffer(m);
		std::string prefix = userInput.prefix; // loads the first map
		std::string firstFile = prefix + "/.map." + std::to_string(m) + ".0.tmp.bin";
		maps[m] = new ParallelMap(0, KeyHasher(seqBuf[m].data, k), KeyEqualTo(seqBuf[m].data, k));
		maps[m]->reserve((uint64_t)(mapSizeCounts[m]*2.5));
		//std::cout<<maps[m]->bucket_count()<<std::endl;
		maps32[m] = new ParallelMap32(0, KeyHasher(seqBuf[m].data, k), KeyEqualTo(seqBuf[m].data, k));
		
		uint8_t fileNum = 0;
		while (fileExists(prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) + ".tmp.bin")) { // for additional map loads the map and merges it
			std::string nextFile = prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) + ".tmp.bin"; // loads the next map
			ParallelMap *nextMap = new ParallelMap(0, KeyHasher(seqBuf[m].data, k), KeyEqualTo(seqBuf[m].data, k));
			phmap::BinaryInputArchive ar_in(nextFile.c_str());
			nextMap->phmap_load(ar_in);
			if (!userInput.keepTmp)
				remove(nextFile.c_str());
			maps[m]->insert(nextMap->begin(), nextMap->end());
			delete nextMap;
			maps32[m]->insert(tmpMaps32[m][fileNum]->begin(), tmpMaps32[m][fileNum]->end());
			delete tmpMaps32[m][fileNum];
			++fileNum;
		}
		summary(m);
		delete seqBuf[m].data;
		dumpMap(userInput.prefix, m);
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
	for(std::thread& thread : writeThreads)
		thread.join();
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
bool Kmap<DERIVED, INPUT, KEY, TYPE1, TYPE2>::dumpTmpMap(std::string prefix, uint8_t m, ParallelMap *map, uint8_t fileNum) {
	
	prefix.append("/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin");
	phmap::BinaryOutputArchive ar_out(prefix.c_str());
	map->phmap_dump(ar_out);
	delete map;
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
		lg.verbose("Memory in use/allocated/total: " + std::to_string(get_mem_inuse(3)) + "/" + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
	
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
	maps[m] = new ParallelMap(0, KeyHasher(seqBuf[m].data, k), KeyEqualTo(seqBuf[m].data, k));
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
	
	maps32[m] = new ParallelMap32(0, KeyHasher(seqBuf[m].data, k), KeyEqualTo(seqBuf[m].data, k));
	
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
	
	if (userInput.keepTmp)
		return;
	
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
