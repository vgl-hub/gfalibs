#ifndef BED_H
#define BED_H

#include <vector>
#include <string>

#include "parallel-hashmap/phmap.h"

class BedCoordinates { // generic representation of bed coordinates
private:
    
    phmap::flat_hash_map<std::string,std::vector<std::pair<uint64_t,uint64_t>>> coordinates;
    std::vector<std::string> headers;
    
public:
    
    void pushCoordinates(std::string h, uint64_t s = 0, uint64_t e = 0) { // reading coordinates
        
        std::pair<uint64_t,uint64_t> newCoordinates(s,e);
        auto got = coordinates.find(h); // insert or find this header in the hash table
        
        if (got == coordinates.end()){
            headers.push_back(h);
            std::vector<std::pair<uint64_t,uint64_t>> newCoordinatesVector = {newCoordinates};
            std::pair<std::string,std::vector<std::pair<uint64_t,uint64_t>>> newCoordinatesVectorPair = std::make_pair(h,newCoordinatesVector);
            coordinates.insert(newCoordinatesVectorPair);
        }else{
            got->second.push_back(newCoordinates);
        }
    }
    
    inline bool empty() { // check if no coordinates are present
        return (coordinates.size()==0) ? true : false; // check if no coordinates are present
    }
    
    inline uint32_t size() { // check if no coordinates are present
        return coordinates.size(); // check if no coordinates are present
    }
    
    inline std::vector<std::string> getHeaders() { // get all the headers
        return headers;
    }
    
    inline bool isPresent(std::string h) { // get a specific header
        auto got = coordinates.find(h); // insert or find this kmer in the hash table
        
        if (got == coordinates.end())
            return false;
        else
            return true;
    }
    
    inline phmap::flat_hash_map<std::string,std::vector<std::pair<uint64_t,uint64_t>>>& getCoordinates() { // get a specific start coordinate
        return coordinates; // insert or find this kmer in the hash table
    }
    
};

#endif /* BED_H */
