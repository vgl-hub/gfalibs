#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

struct membuf : std::streambuf
{
    void set(char* gbeg, char* gnext, char* gend) {
        this->setg(gbeg, gnext, gend);
    }

};

class StreamObj {
    
    std::streambuf* buffer;
    std::shared_ptr<std::istream> stream;
    std::ifstream ifs;
    zstream::igzstream zfin, zin;
    bool file = false, gzip = false, decompress = true;
    unsigned int bufSize = 1000000;
    char* bufContent = NULL, *contents = NULL;
    std::condition_variable mutexCondition;
    membuf sbuf; // streambuffer to hold gzip buffer
    
public:
    
    StreamObj() :
    zfin(ifs), zin(std::cin) {}
    
    ~StreamObj(){this->closeStream(); delete[] bufContent;}
    
    bool isGzip(std::streambuf* buffer);
    
    void decompressBuf(gzFile buffer);

    void readBuf();
    
    std::shared_ptr<std::istream> openStream(UserInput &userInput, char type, unsigned int* file = NULL);
    
    void closeStream();
    
    std::string type();
    
};

#endif /* STREAM_OBJ_H */
