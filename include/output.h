//
//  gfastats-output.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef OUTPUT_H
#define OUTPUT_H

#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"

#include "zlib.h"
#include "zstream/zstream_common.hpp"
#include "zstream/ozstream.hpp"
#include "zstream/ozstream_impl.hpp"

// struct
struct OutputStream {
    
    std::string file;
    std::unique_ptr<std::ostream> stream; // here we create a smart pointer to handle any kind of output stream
    bool gzip = false, outFile = false; // variables to handle output type
    
    std::ofstream ofs; // this stream outputs to file
    zstream::ogzstream zfout; // this stream outputs gzip compressed to file
    zstream::ogzstream zout; // this stream outputs gzip compressed to stdout
    
    OutputStream(std::string file);
    ~OutputStream();
};

//classes
class Report {

unsigned int counter = 0;
    
public:
    bool segmentReport(InSequences &inSequences, int &outSequence_flag);
    
    bool pathReport(InSequences &inSequences);
    
    void writeToStream(InSequences &inSequences, std::string file, UserInput &userInput);
    
    bool outCoord(InSequences &inSequences, char bedOutType, bool sizeOnly = false);
    
    bool reportStats(InSequences &inSequences, uint64_t gSize, int outBubbles_flag);
    
    bool nstarReport(InSequences &inSequences, uint64_t gSize);
    
};

#endif /* OUTPUT_H */
