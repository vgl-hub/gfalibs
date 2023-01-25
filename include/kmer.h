#ifndef KMER
#define KMER

#include <parallel_hashmap/phmap.h>
#include "parallel_hashmap/phmap_dump.h"

struct buf64 {
    uint64_t pos = 0, size = 100000;
    uint64_t *seq = new uint64_t[size];
};

template<class INPUT, typename VALUE>
class Kmap {
    
    std::vector<Log> logs;
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
    
    unsigned int batchSize = 10000;

    uint8_t k;
    
    uint64_t totKmers = 0, totKmersUnique = 0, totKmersDistinct = 0;
    
    const uint64_t mapCount = k < 28 ? pow(4,k/4) : pow(4,6);
    
    const uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount;
    
    uint64_t* pows = new uint64_t[k];

    std::vector<buf64*> buffers;
    
    phmap::flat_hash_map<uint64_t, VALUE>* map = new phmap::flat_hash_map<uint64_t, VALUE>[mapCount];
    
    phmap::flat_hash_map<uint64_t, uint64_t> histogram1, histogram2;
    
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
    
public:
    
    Kmap(uint8_t k) : k(k) {
        
        for(uint8_t p = 0; p<k; ++p)
            pows[p] = (uint64_t) pow(4,p);
        
    };
    
    ~Kmap(){
        
        delete[] map;
        delete[] pows;
        
    }
    
    void convert(INPUT& userInput);
    
    void load(INPUT& userInput);
    
    bool traverseInReads(Sequences* readBatch);
    
    void appendReads(Sequences* readBatch);
    
    inline uint64_t hash(uint8_t* string);
    
    void hashSequences(Sequences* readBatch);
    
    void hashSegments(std::vector<InSegment*>* segments);
    
    void count();
    
    bool countBuff(uint16_t m);
    
    bool histogram(phmap::flat_hash_map<uint64_t, VALUE>& map);
    
    void resizeBuff(buf64* buff);
    
    void printHist(std::unique_ptr<std::ostream>& ostream);
    
    void report(INPUT& userInput);
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);

};

template<class INPUT, typename VALUE>
void Kmap<INPUT, VALUE>::convert(INPUT& userInput) {
    
    InSequences inSequences;
    
    if (!userInput.iSeqFileArg.empty() || userInput.pipeType == 'f') {
        
        StreamObj streamObj;
        
        stream = streamObj.openStream(userInput, 'f');
        
        if (stream) {
            
            switch (stream->peek()) {
                    
                case '>': {
                    
                    stream->get();
                    
                    while (getline(*stream, newLine)) {
                        
                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }
                        
                        std::string* inSequence = new std::string;
                        
                        getline(*stream, *inSequence, '>');
                        
                        lg.verbose("Individual fasta sequence read");
                        
                        Sequence* sequence = new Sequence{seqHeader, seqComment, inSequence, NULL};
                        
                        sequence->seqPos = seqPos; // remember the order
                        
                        inSequences.appendSequence(sequence);
                        
                        seqPos++;
                        
                    }
                    
                    jobWait(threadPool);
                    
                    std::vector<Log> logs = inSequences.getLogs();
                    
                    //consolidate log
                    for (auto it = logs.begin(); it != logs.end(); it++) {
                        
                        it->print();
                        logs.erase(it--);
                        
                    }
                    
                    hashSegments(inSequences.getInSegments());
                    
                    break;
                    
                }
                    
                case '@': {
                    
                    Sequences* readBatch = new Sequences;

                    while (getline(*stream, newLine)) { // file input

                        newLine.erase(0, 1);

                        h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        seqHeader = h;
                        
                        if (c != NULL) {
                            
                            seqComment = std::string(c);
                            
                        }

                        std::string* inSequence = new std::string;
                        getline(*stream, *inSequence);

                        getline(*stream, newLine);
                        
                        ignore(*stream, '\n');

                        readBatch->sequences.push_back(new Sequence {seqHeader, seqComment, inSequence});
                        seqPos++;

                        if (seqPos % batchSize == 0) {

                            readBatch->batchN = seqPos/batchSize;
                            
                            lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                            appendReads(readBatch);

                            readBatch = new Sequences;

                        }

                        lg.verbose("Individual fastq sequence read: " + seqHeader);

                    }
                    
                    readBatch->batchN = seqPos/batchSize + 1;
                        
                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                    appendReads(readBatch);
                    
                    jobWait(threadPool);
                    
                    std::vector<Log> logs = inSequences.getLogs();
                    
                    //consolidate log
                    for (auto it = logs.begin(); it != logs.end(); it++) {
                        
                        it->print();
                        logs.erase(it--);
                        
                    }

                    break;

                }
                    
            }
            
        }
        
    }
    
}

template<class INPUT, typename VALUE>
void Kmap<INPUT, VALUE>::load(INPUT& userInput){
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return loadMap(userInput.iSeqFileArg, m); });
    
    jobWait(threadPool);
    
    lg.verbose("Regenerate histogram");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return histogram(map[m]); });
    
    jobWait(threadPool);
    
}

template<class INPUT, typename VALUE>
bool Kmap<INPUT, VALUE>::loadMap(std::string prefix, uint16_t m) { // loading prototype
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    map[m].phmap_load(ar_in);
    
    return true;

}

template<class INPUT, typename VALUE>
bool Kmap<INPUT, VALUE>::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    map[m].phmap_dump(ar_out);
    
    return true;
    
}

template<class INPUT, typename VALUE>
void Kmap<INPUT, VALUE>::report(INPUT& userInput) {
    
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

template<class INPUT, typename VALUE>
bool Kmap<INPUT, VALUE>::traverseInReads(Sequences* readBatch) { // traverse the read

    hashSequences(readBatch);
    
    delete readBatch;
    
    return true;
    
}

template<class INPUT, typename VALUE>
void Kmap<INPUT, VALUE>::appendReads(Sequences* readBatch) { // read a collection of reads
    
    threadPool.queueJob([=]{ return traverseInReads(readBatch); });
    
    std::unique_lock<std::mutex> lck(mtx);
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        
    }
    
}

template<class INPUT, typename VALUE>
inline uint64_t Kmap<INPUT, VALUE>::hash(uint8_t *kmer) {
    
    uint64_t fw = 0, rv = 0;
    
    for(uint8_t c = 0; c<k; ++c)
        fw += *kmer++ * pows[c];
    
    --kmer;
    
    for(uint8_t c = 0; c<k; ++c)
        rv += (3-(*kmer--)) * pows[c];
    
    return fw < rv ? fw : rv;
}

template<class INPUT, typename VALUE>
bool Kmap<INPUT, VALUE>::countBuff(uint16_t m) {

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    buf64* thisBuf;
    
    phmap::flat_hash_map<uint64_t, VALUE>* thisMap;
    
    for(buf64* buf : buffers) {
        
        thisBuf = &buf[m];
        
        thisMap = &map[m];
        
        uint64_t len = thisBuf->pos;
        
        for (uint64_t c = 0; c<len; ++c)
            ++(*thisMap)[thisBuf->seq[c]];
        
        delete[] thisBuf->seq;
        
    }
    
    return true;

}

template<class INPUT, typename VALUE>
bool Kmap<INPUT, VALUE>::histogram(phmap::flat_hash_map<uint64_t, VALUE>& map) {
    
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
        
        histogram1[pair.first] += pair.second;
        
        totKmers += pair.first * pair.second;
        
    }
    
    return true;

}

template<class INPUT, typename VALUE>
void Kmap<INPUT, VALUE>::printHist(std::unique_ptr<std::ostream>& ostream) {
    
    std::vector<std::pair<uint64_t, uint64_t>> table(histogram1.begin(), histogram1.end());
    std::sort(table.begin(), table.end());
    
    for (auto pair : table)
        *ostream<<pair.first<<"\t"<<pair.second<<"\n";

}

template<class INPUT, typename VALUE>
void Kmap<INPUT, VALUE>::hashSequences(Sequences* readBatch) {
    
    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    buf64* buf = new buf64[mapCount];
    
    for (Sequence* sequence : readBatch->sequences) {
        
        uint64_t len = sequence->sequence->size(), Kmap = len-k+1;
        
        if (len<k)
            continue;
        
        unsigned char* first = (unsigned char*)sequence->sequence->c_str();
        
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        uint64_t value, i, newSize;
        buf64* b;
        uint64_t* bufNew;
        
        for (uint64_t c = 0; c<Kmap; ++c){
            
            value = hash(str+c);
            
            i = value / moduloMap;
            
            b = &buf[i];
            
            if (b->pos == b->size) {
                
                newSize = b->size * 2;
                bufNew = new uint64_t[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            b->seq[b->pos++] = value;
                        
        }
        
        delete[] str;
        
        threadLog.add("Processed sequence: " + sequence->header);
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
    logs.push_back(threadLog);
    
}

template<class INPUT, typename VALUE>
void Kmap<INPUT, VALUE>::hashSegments(std::vector<InSegment*>* segments) {
    
    buf64* buf = new buf64[mapCount];
    
    for (InSegment* segment : *segments) {
        
        uint64_t len = segment->getSegmentLen(), Kmap = len-k+1;
        
        if (len<k)
            continue;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        uint64_t value, i, newSize;
        buf64* b;
        uint64_t* bufNew;
        
        for (uint64_t c = 0; c<Kmap; ++c){
            
            value = hash(str+c);
            
            i = value / moduloMap;
            
            b = &buf[i];
            
            if (b->pos == b->size) {
                
                newSize = b->size * 2;
                bufNew = new uint64_t[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            b->seq[b->pos++] = value;
                        
        }
        
        delete[] str;
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
}

template<class INPUT, typename VALUE>
void Kmap<INPUT, VALUE>::count() {
    
    lg.verbose("Counting with " + std::to_string(mapCount) + " maps");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return countBuff(m); });
    
    jobWait(threadPool);
    
    lg.verbose("Generate histogram");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return histogram(map[m]); });
    
    jobWait(threadPool);
    
}


#endif //KMER
