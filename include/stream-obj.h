#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

#include "zlib.h"

struct membuf : std::streambuf {
    
    static const unsigned int bufSize = 1000000;
    unsigned int size = bufSize;
    char bufContent1[bufSize], bufContent2[bufSize];
    char* bufContent = bufContent1;
    gzFile fi;
    bool decompressed1 = false, decompressed2 = false, start = false, eof = false, uflowDone1 = true, uflowDone2 = true;
    std::mutex semMtx;
    std::condition_variable semaphore;
    
    void openFile(std::string file);
    
    void read();
    
    int uflow();
    
    bool decompressBuf();
    
};

class StreamObj {
    
    std::streambuf* buffer;
    std::shared_ptr<std::istream> stream;
    std::ifstream ifs;
//    zstream::igzstream zfin, zin;
    bool file = false, gzip = false;
//    std::istringstream strm;
    char* content = new char[10000];
    std::condition_variable mutexCondition;
//    bool decompress = true, done = false;
    
public:
    
//    StreamObj() :
//    zfin(ifs), zin(std::cin) {}
    
    ~StreamObj(){this->closeStream(); delete[] content;}
    
    bool isGzip(std::streambuf* buffer);
    
//    void decompressBuf(std::streambuf* buffer);
//
//    void readBuf(std::streambuf* buffer);
    
    std::shared_ptr<std::istream> openStream(UserInput &userInput, char type, unsigned int* file = NULL);
    
    void closeStream();
    
    std::string type();
    
    std::shared_ptr<std::istream> returnStream();
    
};

#endif /* STREAM_OBJ_H */
