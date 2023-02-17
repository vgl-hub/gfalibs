#ifndef FASTX_H
#define FASTX_H

template<typename OBJECT>
bool loadSequences(UserInput userInput, OBJECT* object, char type, unsigned int* fileNum) { // load from FASTA/FASTQ to templated object
    
    // stream read variables
    char* c;
    std::string firstLine, newLine, seqHeader, seqComment, line, bedHeader;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    unsigned int batchSize = 10000; // number of sequences processed by a thread
        
    if (!userInput.iSeqFileArg.empty() || userInput.pipeType == type) {
                
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
                        Sequence* sequence = new Sequence{seqHeader, seqComment, inSequence, NULL};
                        sequence->seqPos = seqPos; // remember the order
                        object->appendSequence(sequence);
                        seqPos++;
                        
                    }
                    
                    jobWait(threadPool);
                    
                    break;
                    
                }
                    
                case '@': {
                    
                    Sequences* readBatch = new Sequences;

                    while (getline(*stream, newLine)) { // file input

                        newLine.erase(0, 1);
                        seqHeader = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                        c = strtok(NULL,""); //read comment
                        
                        if (c != NULL)
                            seqComment = std::string(c);

                        std::string* inSequence = new std::string;
                        getline(*stream, *inSequence);
                        getline(*stream, newLine);
                        ignore(*stream, '\n');
                        readBatch->sequences.push_back(new Sequence {seqHeader, seqComment, inSequence});
                        seqPos++;

                        if (seqPos % batchSize == 0) {

                            readBatch->batchN = seqPos/batchSize;
                            lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                            threadPool.queueJob([=]{ return object->traverseInReads(readBatch); });
                            std::unique_lock<std::mutex> lck(mtx);
                            for (auto it = object->logs.begin(); it != object->logs.end(); it++) {
                             
                                it->print();
                                object->logs.erase(it--);
                                
                            }
                            
                            readBatch = new Sequences;

                        }

                        lg.verbose("Individual fastq sequence read: " + seqHeader);

                    }
                    
                    readBatch->batchN = seqPos/batchSize + 1;
                    lg.verbose("Processing batch N: " + std::to_string(readBatch->batchN));

                    threadPool.queueJob([=]{ return object->traverseInReads(readBatch); });
                    std::unique_lock<std::mutex> lck(mtx);
                    for (auto it = object->logs.begin(); it != object->logs.end(); it++) {
                     
                        it->print();
                        object->logs.erase(it--);
                        
                    }

                    break;

                }
                    
            }
            
            //consolidate log
            jobWait(threadPool);
            for (auto it = object->logs.begin(); it != object->logs.end(); it++) {
                
                it->print();
                object->logs.erase(it--);
                
            }
            
        }
        
    }
    
    return true;

}

#endif /* FASTX_H */
