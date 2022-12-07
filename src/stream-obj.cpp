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


    
void membuf::openFile(std::string file) {
    
//        std::cout<<"file open"<<std::endl;
    
    fi = gzopen(file.c_str(), "rb");
    threadPool.queueJob([=]{ return decompressBuf(); });
    
    read();
    
}

void membuf::read() {
    
//        std::cout<<"read"<<std::endl;
    
    std::unique_lock<std::mutex> lck(mtx);
    
    mutexCondition.wait(lck, [this] {
        return !decompress;
    });
    
    decompress = false;
    mutexCondition.notify_one();
    
}

char membuf::snextc() {
    
    std::cout<<"snextc"<<std::endl;
    
    if ( sbumpc() == EOF ) return EOF;
    else return sgetc();
    
}

int membuf::sbumpc() {
    
    gbump(1);
    
    if ( (!gptr()) || (gptr()==egptr()) ) {
        
        std::cout<<"resetting buffer"<<std::endl;
        
        decompress = true;
        mutexCondition.notify_one();
        
        std::unique_lock<std::mutex> lck(mtx);
        
        mutexCondition.wait(lck, [this] {
            return !decompress || eof;
        });
        
    }
    
    return gptr()[-1];
    
}

bool membuf::decompressBuf() {

    std::unique_lock<std::mutex> lck(mtx);
    
    unsigned int size = bufSize;
    
    while(size==bufSize) {
        
        mutexCondition.wait(lck, [this] {
//                std::cout<<"decompression thread is waiting"<<std::endl;
            return decompress;
        });
        
        size = gzread(fi, bufContent, sizeof(char)*bufSize);
        
        set(bufContent, bufContent, bufContent + sizeof(bufContent));
        
        decompress = false;
        
        mutexCondition.notify_one();
        
    };
    
    eof = true;
    
    gzclose(fi);
    
    return eof;

}

membuf* memstream::rdbuf() {return assBuf;}

bool getline(memstream& in, std::string& newLine) {

    char c = in.rdbuf()->sgetc();
    
    newLine.clear();

     do{

        newLine.push_back(c);

     }while ((c = in.rdbuf()->membuf::snextc()) != '\n' && c != EOF);

    return c != EOF ? true : false;

}

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
            
            membuf sbuf;
            
            sbuf.openFile(userInput.file(type));
            
            memstream in(&sbuf);
            
            return std::make_shared<memstream>(&sbuf);
                

        }

    }else{
        
        buffer = std::cin.rdbuf();
        
        gzip = isGzip(buffer);

        if (gzip) {
            
            //TBD

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
