#ifndef BED_H
#define BED_H

#include <vector>
#include <string>

class BedCoordinates { // generic representation of bed coordinates
private:
    std::vector<std::string> seqHeaders;
    std::vector<unsigned int> cBegin;
    std::vector<unsigned int> cEnd;
    
public:
    
    void pushCoordinates(std::string h, unsigned int b = 0, unsigned int e = 0) { // reading coordinates
        seqHeaders.push_back(h);
        cBegin.push_back(b);
        cEnd.push_back(e);
    }
    
    inline bool empty() { // check if no coordinates are present
        return (seqHeaders.size()==0) ? true : false; // check if no coordinates are present
    }
    
    inline unsigned int size(){
        return seqHeaders.size(); // check if no coordinates are present
    }
    
    inline std::vector<std::string> getSeqHeaders() { // get all the headers
        return seqHeaders;
    }
    
    inline std::string getSeqHeader(unsigned int pos) { // get a specific header
        return seqHeaders[pos];
    }
    
    inline unsigned int getcBegin(unsigned int pos) { // get a specific start coordinate
        return cBegin[pos];
    }
    
    inline unsigned int getcEnd(unsigned int pos) { // get a specific end coordinate
        return cEnd[pos];
    }
    
};

#endif /* BED_H */
