#ifndef FASTX_H
#define FASTX_H

#include "stream-obj.h"

template<typename OBJECT>
bool loadSequences(UserInput userInput, OBJECT& object, char type, uint16_t fileNum) { // load from FASTA/FASTQ to templated object
    
    // stream read variables
    char* c;
    std::string firstLine, newLine, seqHeader, seqComment, line, bedHeader;
    uint32_t seqPos = 0; // to keep track of the original sequence order
    uint32_t batchSize = 10000; // number of sequences processed by a thread
                
    //stream objects
    StreamObj streamObj;
    std::shared_ptr<std::istream> stream;
    stream = streamObj.openStream(userInput, type, fileNum);
    
    if (stream) {
        
        switch (stream->peek()) {
                
            case '>': {
                
                stream->get();
                
                while (getline(*stream, newLine)) {
                    
                    seqHeader = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                    c = strtok(NULL,""); //read comment
                    
                    if (c != NULL)
                        seqComment = std::string(c);
                    
                    std::string* inSequence = new std::string;
                    getline(*stream, *inSequence, '>');
                    lg.verbose("Individual fasta sequence read");
                    Sequence* sequence = new Sequence{seqHeader, seqComment, inSequence, NULL, seqPos++};
                    object.appendSequence(sequence);
                }
                break;
            }
                
            case '@': {
                
                Sequences* readBatch = new Sequences;

                while (getline(*stream, newLine)) { // file input
                    
                    newLine.erase(0, 1);
                    seqHeader = std::string(strtok(strdup(newLine.c_str())," ")); // process header line
                    c = strtok(NULL,""); // read comment
                    
                    if (c != NULL)
                        seqComment = std::string(c);

                    std::string* inSequence = new std::string;
                    getline(*stream, *inSequence);
                    getline(*stream, newLine);
                    ignore(*stream, '\n');
                    readBatch->sequences.push_back(new Sequence {seqHeader, seqComment, inSequence});
                    ++seqPos;

//                    lg.verbose("Individual fastq sequence read: " + seqHeader);
                    
                    if (seqPos % batchSize == 0) {
                        readBatch->batchN = seqPos/batchSize;
                        object.traverseInReads(readBatch);
                        object.consolidate();
                        readBatch = new Sequences;
                    }
                }
                readBatch->batchN = seqPos/batchSize+1;
                object.traverseInReads(readBatch);
                break;
            }
        }
    }
    return true;
}

template<typename OBJECT>
bool loadKmers(UserInput userInput, OBJECT& object, char type, uint16_t fileNum) { // load from FASTA/FASTQ to templated object, faster when we only need to retain kmers not the original reads
    
    // stream read variables
    int batchSize = 10000; // number of bases processed
                
    //stream objects
    StreamObj streamObj;
    std::shared_ptr<std::istream> stream;
    stream = streamObj.openStream(userInput, type, fileNum);
    
    if (stream) {
        
        switch (stream->peek()) {
                
            case '>': {
                
                
                break;
                
            }
                
            case '@': {
                
                std::string* readBatch;
                
                while (*stream) { // file input
                    allocMemory(batchSize * sizeof(char));
                    readBatch = new std::string;
                    getKmers(*stream, *readBatch, batchSize);
                    object.traverseInReads(readBatch);
                    object.consolidate();
                }
                break;
            }
        }
    }
    return true;
}

#endif /* FASTX_H */
