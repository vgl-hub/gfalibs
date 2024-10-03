#include <cstdint>
#include <stdlib.h>
#include <fstream>
#include <vector>

#include "bed.h"
#include "struct.h"
#include "zlib.h"
#include "global.h"
#include "log.h"
#include "threadpool.h"

#include "stream-obj.h"

void membuf::openFile(std::string file) {
    
//    std::cout<<"file open: "<<file<<std::endl;
    
    fi = gzopen(file.c_str(), "rb");
    
    decompressor = new std::thread(&membuf::decompressBuf, this);
    
    wait();
    
}

void membuf::wait() {
    
//    std::cout<<"R:wait"<<std::endl;
    
    {
        
        std::unique_lock<std::mutex> lck(semMtx);
        
        semaphore.wait(lck, [this] {
            return start;
        });
        
    }
    
//    std::cout<<"R:read"<<std::endl;
    
}

int membuf::uflow() {
    
//    std::cout<<"R:resetting buffer"<<std::endl;
    
    {
        
//        std::cout<<"R:waiting for decompressed buffer"<<std::endl;
        
        std::unique_lock<std::mutex> lck(semMtx);
        
        semaphore.wait(lck, [this] {
            return (decompressed1 && whichBuf) || (decompressed2 && !whichBuf) || eof;
        });
    
        if(decompressed1 && whichBuf) {
            
//            std::cout<<"R:setting internal buffer to buffer 1"<<std::endl;
            
            setg(bufContent1, bufContent1, bufContent1 + sizeof(char)*size1);
            
            decompressed1 = false;
            whichBuf = 0;
            
        }
        
        if (decompressed2 && !whichBuf){
            
//            std::cout<<"R:setting internal buffer to buffer 2"<<std::endl;
            
            setg(bufContent2, bufContent2, bufContent2 + sizeof(char)*size2);
            
            decompressed2 = false;
            whichBuf = 1;

        }
        
    }
    
    semaphore.notify_one();
    
    if (sgetc() == EOF) {
        if (decompressor != NULL && decompressor->joinable()) {
            decompressor->join();
            delete decompressor;
            decompressor = NULL;
        }
        return EOF;
    }
    gbump(1);
    return gptr()[-1];
    
}

bool membuf::decompressBuf() {
    
    *size = gzread(fi, bufContent, sizeof(char)*bufSize);
    
    setg(bufContent, bufContent, bufContent + sizeof(char)**size);
    
    start = true;
    
    semaphore.notify_one();
    
//    std::cout<<"D:extracted bases: "<<*size<<std::endl;
    
    while(*size==bufSize) {
         
        bufContent = (bufContent == bufContent1) ? bufContent2 : bufContent1;
        
        size = (bufContent == bufContent1) ? &size1 : &size2;
        
//        std::cout<<"D:buffer swapped"<<std::endl;
        
        *size = gzread(fi, bufContent, sizeof(char)*bufSize);
        
        {
            
            std::unique_lock<std::mutex> lck(semMtx);
            
            (bufContent == bufContent1) ? decompressed1 = true : decompressed2 = true;
            
        }
        
        semaphore.notify_one();
        
//        std::cout<<"D:extracted bases: "<<*size<<std::endl;
        
        {
            
//            std::cout<<"D:waiting for buffer being read"<<std::endl;
            
            std::unique_lock<std::mutex> lck(semMtx);
            
            semaphore.wait(lck, [this] {
                
                return (bufContent == bufContent1) ? !whichBuf : whichBuf;
                
            });
            
        }
        
    };
    
    eof = true;
    
    semaphore.notify_one();
    
    gzclose(fi);
    
//    std::cout<<"D:decompression completed"<<std::endl;
    
    return eof;

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

std::shared_ptr<std::istream> StreamObj::openStream(UserInput& userInput, char type, uint16_t fileNum) {
    
    file = userInput.pipeType != type ? true : false;
    
    if (file) { // input is from file
        
        if (fileNum == 0) {
        
            ifs.open(userInput.file(type)); // this stream takes input from a plain file
            
        }else{
            
            ifs.open(userInput.file(type, fileNum));
            
        }

        buffer = ifs.rdbuf();
        
        gzip = isGzip(buffer);

        if (gzip) {
            
            sbuf.openFile(userInput.file(type, fileNum));
            
            buffer = &sbuf;

        }

    }else{
        
        buffer = std::cin.rdbuf();
        
        gzip = isGzip(buffer);

        if (gzip) {
            
            //TBD

        }
        
    }
    
    lg.verbose("Created stream object from input.\nStream type (" + this->type() + ").\nStreaming started.");
    
    return std::make_shared<std::istream>(buffer);
    
}

void StreamObj::closeStream() {

    if (gzip) {
    }
        
    ifs.close();
        
    lg.verbose("File stream closed");

}
