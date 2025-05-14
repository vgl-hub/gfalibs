#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

#include <fstream>
#include "zlib.h"

class membuf : public std::streambuf {
    
#ifdef _WIN32
    static const unsigned int bufSize = 500000;
#else
    static const unsigned int bufSize = 1000000;
#endif
    unsigned int size1 = bufSize, size2 = 0;
    unsigned int* size = &size1;
    char bufContent1[bufSize], bufContent2[bufSize];
    char* bufContent = bufContent1;
    gzFile fi;
    bool decompressed1 = false, decompressed2 = false, start = false, eof = false, whichBuf = 0;
    std::mutex semMtx;
    std::condition_variable semaphore;
    std::thread *decompressor;
    
public:
    
    void openFile(std::string file);
    
    void wait();
    
    int uflow();
    
    bool decompressBuf();
	
	void close();
    
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
    
    std::shared_ptr<std::istream> openStream(UserInput &userInput, char type, uint16_t file = 0);
    
    void closeStream();
    
    std::string type();
    
};

#endif /* STREAM_OBJ_H */
