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
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "global.h"
#include "log.h"

#include "threadpool.h"
#include "stream-obj.h"

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

//void StreamObj::decompressBuf(std::streambuf* buffer) {
//    
//    int chars = 0;
//    
//    std::unique_lock<std::mutex> lck(mtx, std::defer_lock);
//    
//    while (true) {
//        
//        if (decompress) {
//        
//            lck.lock();
//                                    
//            chars = buffer->sgetn(content, 1024);
//            std::cout<<"decompress: "<<chars<<" "<<std::endl;
//
//            strm.str(content);
//            
//            decompress = false;
//            
//            if (chars < 1024) {
//                done = true;
//                return;
//            }
//            
//            lck.unlock();
//            
//        }
//    
//    }
//    
//}
//
//void StreamObj::readBuf(std::streambuf* buffer) {
//    
//    int chars = 0;
//    
//    std::unique_lock<std::mutex> lck(mtx, std::defer_lock);
//    
//    while (true) {
//        
//        if (!decompress) {
//        
//            lck.lock();
//            
//            std::string line;
//
//            getline(strm,line);
//            
//            std::cout<<line<<std::endl;
//
//            if (!strm.eof())
//                decompress = true;
//            else
//                return;
//            
//            if (done)
//                decompress = false;
//            
//            lck.unlock();
//            
//        }
//        
//    }
//    
//}

std::shared_ptr<std::istream> StreamObj::openStream(UserInput& userInput, char type, unsigned int* fileNum) {
    

            
            zfin.open();
            
//            threadPool.queueJob([=]{ return decompressBuf(zfin.rdbuf()); });
            
//            buffer = strm.rdbuf();
            
            buffer = zfin.rdbuf();
            
            do {
              char ch = buffer->sbumpc();
              std::cout << ch;
              if(buffer->snextc() == EOF)
                std::cout << "EOF present\n";
                
            } while (true);
            
            std::cout << "im done\n";
            
            exit(1);


    
    return std::make_shared<std::istream>(buffer);
    
}

std::shared_ptr<std::istream> StreamObj::returnStream() {
    
    return stream;
    
}

void StreamObj::closeStream() {

    if (gzip) {

        zfin.read_footer();

        if (zfin.check_crc()) {

            lg.verbose("Crc check successful");

        }else{

            lg.verbose("Warning: crc check unsuccessful. Check input file");

        }

    }
        
    ifs.close();
        
    lg.verbose("File stream closed");

}
