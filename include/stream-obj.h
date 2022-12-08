#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

#include "zlib.h"

struct membuf : std::streambuf {
    
    unsigned int bufSize = 1000000;
    char bufContent[1000000];
    gzFile fi;
    bool decompress = true, eof = false;
    std::mutex semMtx;
    std::condition_variable mutexCondition;
    
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
