#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

#include <fstream>
#include "zlib.h"

class membuf : public std::streambuf {
    
#ifdef _WIN32
    static const unsigned int bufSize = 500000;
#else
    static const unsigned int bufSize = 1048576;
#endif

    unsigned int size1 = bufSize, size2 = 0;
    unsigned int* size = &size1;
	std::unique_ptr<char[]> bufContent1{ new char[bufSize] };
	std::unique_ptr<char[]> bufContent2{ new char[bufSize] };
    char* bufContent = nullptr;
    gzFile fi;
    bool decompressed1 = false, decompressed2 = false, start = false, eof = false, whichBuf = 0;
    std::mutex semMtx;
    std::condition_variable semaphore;
    std::thread *decompressor;
	std::streampos pos;
    
public:
	
	membuf() { bufContent = bufContent1.get(); }
    
    void openFile(std::string file);
    
    void wait();
    
    int uflow() override;
    
    bool decompressBuf();
	
	void gzClose();
	
	std::streampos seekoff(std::streamoff off, std::ios_base::seekdir dir, std::ios_base::openmode) override {
		auto pos = gptr();
		if (dir == std::ios_base::cur)
			pos += off;
		else if (dir == std::ios_base::end)
			pos = egptr() + off;
		else if (dir == std::ios_base::beg)
			pos = eback() + off;

		//check bounds
		if (pos < eback())
			return std::streambuf::pos_type(-1);
		else if (pos > egptr())
			return std::streambuf::pos_type(-1);

		setg(eback(), pos, egptr());
		return gptr() - eback();
	}

	std::streampos seekpos(std::streampos sp, std::ios_base::openmode which) override {
		return seekoff(std::streamoff(sp), std::ios_base::beg, which);
	}
    
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
