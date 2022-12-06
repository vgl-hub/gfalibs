#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <cstring>

#include <queue>
#include <thread>
#include <condition_variable>
#include <mutex>

#include "bed.h"
#include "struct.h"

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "global.h"
#include "log.h"

#include "threadpool.h"
#include "stream-obj.h"

std::string StreamObj::type() {
    
    std::string type;
    type =  file ? "file" : "pipe";
    type += "/";
    type += gzip ? "gzip" : "plain";
    
    return type;
    
}

bool StreamObj::isGzip(std::streambuf* buffer) {
    
    unsigned char header[2];
    
    header[0] = buffer->sgetc();
    header[1] = buffer->snextc();
    buffer->sungetc();
    
    return (header[0] == 0x1f && header[1] == 0x8b) ? true : false;
    
}

void StreamObj::decompressBuf(gzFile fi) {
    
    bufContent = new char[bufSize];

    std::unique_lock<std::mutex> lck(mtx);
    
    while(gzread(fi, bufContent, bufSize)) {
        
        mutexCondition.wait(lck, [this] {
            return decompress;
        });
        sbuf.set(bufContent, bufContent, bufContent + bufSize);
        
        decompress = false;
        
        mutexCondition.notify_one();
        
    }
    
    sbuf.set(bufContent, bufContent + bufSize, bufContent + bufSize);

}

void StreamObj::readBuf() {
    
    contents = new char [bufSize];
    
    unsigned int n = bufSize;
    
    std::unique_lock<std::mutex> lck(mtx);
    
    while (n == bufSize) {
        
        mutexCondition.wait(lck, [this] {
            return !decompress;
        });
        
        n = sbuf.sgetn(contents, bufSize);
        
        decompress = true;

        mutexCondition.notify_one();
        
    }
    
}

std::shared_ptr<std::istream> StreamObj::openStream(UserInput& userInput, char type, unsigned int* fileNum) {
    
    file = userInput.pipeType != type ? true : false;
    
    if (file) { // input is from file
        
        ifs.close();
        
        if (fileNum == NULL) {
        
            ifs.open(userInput.file(type)); // this stream takes input from a plain file
            
        }else{
            
            ifs.open(userInput.file(type, fileNum));
            
        }

        buffer = ifs.rdbuf();
        
        gzip = isGzip(buffer);

        if (gzip) {
                        
            gzFile fi = gzopen(userInput.file(type).c_str(), "rb");
            
            threadPool.queueJob([=]{ return decompressBuf(fi); });
            
            readBuf();
            
            buffer = &sbuf;

        }

    }else{
        
        buffer = std::cin.rdbuf();
        
        gzip = isGzip(buffer);

        if (gzip) {
            
            zin.open();
            
            buffer = zin.rdbuf();

        }
        
    }
    
    return std::make_shared<std::istream>(buffer);
    
}

void StreamObj::closeStream() {

    if (gzip) {

//        zfin.read_footer();
//
//        if (zfin.check_crc()) {
//
//            lg.verbose("Crc check successful");
//
//        }else{
//
//            lg.verbose("Warning: crc check unsuccessful. Check input file");
//
//        }

    }
        
    ifs.close();
        
    lg.verbose("File stream closed");

}
