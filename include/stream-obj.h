#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

#include "zlib.h"

class membuf : std::streambuf {
    
    static const unsigned int bufSize = 1000000;
    unsigned int size = bufSize;
    char bufContent1[bufSize], bufContent2[bufSize];
    char* bufContent = bufContent1;
    gzFile fi;
    bool decompressed1 = false, decompressed2 = false, start = false, eof = false, uflowDone1 = true, uflowDone2 = true;
    std::mutex semMtx;
    std::condition_variable semaphore;
    
public:
    
    void openFile(std::string file);
    
    void read();
    
    int uflow();
    
    bool decompressBuf();
    
};

class StreamObj {
    
    std::streambuf* buffer;
    membuf sbuf;
    std::ifstream ifs;
    bool file = false, gzip = false;
    
public:
    
    ~StreamObj(){this->closeStream();}
    
    bool isGzip(std::streambuf* buffer);

    void readBuf();
    
    std::shared_ptr<std::istream> openStream(UserInput &userInput, char type, unsigned int* file = NULL);
    
    void closeStream();
    
    std::string type();
    
};

#endif /* STREAM_OBJ_H */
