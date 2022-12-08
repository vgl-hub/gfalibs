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

#include "global.h"
#include "log.h"

#include "threadpool.h"
#include "stream-obj.h"


void membuf::openFile(std::string file) {
    
    std::cout<<"file open: "<<file<<std::endl;
    
    fi = gzopen(file.c_str(), "rb");
    exit(1);
    threadPool.queueJob([=]{ return decompressBuf(); });
    
    read();
    
}

void membuf::read() {
    
        std::cout<<"read"<<std::endl;
    {
        
        std::unique_lock<std::mutex> lck(semMtx);
        
        mutexCondition.wait(lck, [this] {
            return !decompress;
        });
        
        decompress = false;
        
    }
    
    mutexCondition.notify_one();
    
}

int membuf::uflow() {
    
    std::cout<<"resetting buffer"<<std::endl;
    
    {
            
        std::unique_lock<std::mutex> lck(semMtx);
        
        mutexCondition.wait(lck, [this] {
            return !decompress || eof;
        });
        
        decompress = true;
        
    }
    
    mutexCondition.notify_one();
    
    return gptr()[-1];
    
}

bool membuf::decompressBuf() {
    
    unsigned int size = bufSize;
    
    while(size==bufSize) {
        
        {
            std::unique_lock<std::mutex> lck(semMtx);
        
            mutexCondition.wait(lck, [this] {
                std::cout<<"decompression thread is waiting"<<std::endl;
                return decompress;
            });
            
            size = gzread(fi, bufContent, sizeof(char)*bufSize);
            
            std::cout<<size<<std::endl;
            
            set(bufContent, bufContent, bufContent + sizeof(bufContent) - sizeof(char)*(bufSize-size));
            
            decompress = false;
            
            std::cout<<"buffer replenished"<<std::endl;
        
        }
        
        mutexCondition.notify_one();
        
    };
    
    eof = true;
    
    gzclose(fi);
    
    std::cout<<"decompression completed"<<std::endl;
    
    return eof;

}

membuf* memstream::rdbuf() {return assBuf;}

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
        
        if (fileNum == NULL) {
        
            ifs.open(userInput.file(type)); // this stream takes input from a plain file
            
        }else{
            
            ifs.open(userInput.file(type, fileNum));
            
        }

        buffer = ifs.rdbuf();
        
        gzip = isGzip(buffer);

        if (gzip) {
            
            sbuf.openFile(userInput.file(type));
            
            buffer = &sbuf;
            
            char c;
            
            do{
                
                in.get(c);
                std::cout<<c;
                
            }while (c != EOF);
            
            exit(1);

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
