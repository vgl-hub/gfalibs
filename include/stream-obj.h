#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

struct membuf : std::streambuf
{
    unsigned int bufSize = 1000000;
    char bufContent[1000000];
    gzFile fi;
    bool decompress = true, eof = false;
    std::condition_variable mutexCondition;
    
    void set(char* gbeg, char* gnext, char* gend) {
        this->setg(gbeg, gnext, gend);
    }
    
    void openFile(std::string file) {
        
        std::cout<<"file open"<<std::endl;
        
        fi = gzopen(file.c_str(), "rb");
        threadPool.queueJob([=]{ return decompressBuf(); });
        
        read();
        
    }
    
    void read() {
        
        std::cout<<"read"<<std::endl;
        
        std::unique_lock<std::mutex> lck(mtx);
        
        mutexCondition.wait(lck, [this] {
            return !decompress;
        });
        
        decompress = false;
        mutexCondition.notify_one();
        
        lck.unlock();
        
        
    }
    
    char snextc() {
        
        if ( sbumpc() == EOF ) return EOF;
        else return sgetc();
        
    }
    
    int sbumpc() {
        
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
    
    void decompressBuf() {

        std::unique_lock<std::mutex> lck(mtx);
        
        while(gzread(fi, bufContent, sizeof(char)*bufSize)) {
            
            mutexCondition.wait(lck, [this] {
                std::cout<<"decompression thread is waiting"<<std::endl;
                return decompress;
            });
            
            set(bufContent, bufContent, bufContent + sizeof(bufContent));
            
            decompress = false;
            
            mutexCondition.notify_one();
            
        }
        
        eof = true;
        
        gzclose(fi);

    }

};

class StreamObj {
    
    std::streambuf* buffer;
    std::shared_ptr<std::istream> stream;
    std::ifstream ifs;
    zstream::igzstream zfin, zin;
    bool file = false, gzip = false;
    
public:
    
    StreamObj() :
    zfin(ifs), zin(std::cin) {}
    
    ~StreamObj(){this->closeStream();}
    
    bool isGzip(std::streambuf* buffer);
    
    void decompressBuf(gzFile buffer);

    void readBuf();
    
    std::shared_ptr<std::istream> openStream(UserInput &userInput, char type, unsigned int* file = NULL);
    
    void closeStream();
    
    std::string type();
    
};

#endif /* STREAM_OBJ_H */
