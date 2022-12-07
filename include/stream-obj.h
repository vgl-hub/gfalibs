#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

#include "zlib.h"

struct membuf : std::streambuf {
    
    unsigned int bufSize = 1000000;
    char bufContent[1000000];
    gzFile fi;
    bool decompress = true, eof = false;
    std::condition_variable mutexCondition;
    
    void set(char* gbeg, char* gnext, char* gend) {
        this->setg(gbeg, gnext, gend);
    }
    
    void openFile(std::string file);
    
    void read();
    
    char snextc();
    
    int sbumpc();
    
    bool decompressBuf();
    
};

class memstream : public std::istream {

    membuf* assBuf;
    
public:
    
    memstream(membuf* sbuf) : std::istream(sbuf) {
        
        init(sbuf);
        
        assBuf = sbuf;
        
    }
    
    membuf* rdbuf();
    
};

class StreamObj {
    
    std::streambuf* buffer;
    std::shared_ptr<std::istream> stream;
    std::ifstream ifs;
    bool file = false, gzip = false;
    
public:
    
    ~StreamObj(){this->closeStream();}
    
    bool isGzip(std::streambuf* buffer);
    
    void decompressBuf(gzFile buffer);

    void readBuf();
    
    std::shared_ptr<std::istream> openStream(UserInput &userInput, char type, unsigned int* file = NULL);
    
    void closeStream();
    
    std::string type();
    
};

#endif /* STREAM_OBJ_H */
