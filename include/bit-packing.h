#ifndef BIT_PACKING_H
#define BIT_PACKING_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <bitset>

#include "struct.h"

const unsigned char ctoi[256] = { // this converts ACGT>0123 and any other character to 4 in time O(1)
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

const uint8_t itoc[4] = {'A', 'C', 'G', 'T'}; // 0123>ACGT

template <typename TYPE = uint8_t>
struct Buf2bit : Buf<TYPE> {
    
    using Buf<TYPE>::Buf;
    
    inline uint8_t getBase(uint64_t index) const {
        return (this->seq[index/4] >> (6-(index%4)*2)) & 3;
    }
    
    inline std::string substr(uint64_t offset, uint64_t k) {
        
        uint8_t *binaryStrAt = &this->seq[offset/4];
        uint8_t offsetBit = offset % 4;
        
        std::string str;
        uint8_t res = offsetBit + k % 4;

        for (uint64_t i = 0; i < k/4+res/4; ++i) {
            for (int8_t e = 6-offsetBit*2; e >= 0; e = e - 2)
                str.push_back(itoc[(binaryStrAt[i] >> e) & 3]);
            offsetBit = 0;
        }
        for (uint8_t i = 0; i < res % 4; ++i)
            str.push_back(itoc[(binaryStrAt[k/4+res/4] >> (6-i*2)) & 3]);
        
        return str;
    }
    
};

inline uint8_t* bitDecodeToArr(uint8_t *binaryStr, uint64_t offset, uint64_t k) {
    
    uint8_t *str = new uint8_t[k];
    uint8_t *binaryStrAt = &binaryStr[offset/4];
    uint8_t offsetBit = offset % 4;
    
    uint8_t res = offsetBit + k % 4;
    uint64_t pos = 0;
    for (uint64_t i = 0; i < k/4+res/4; ++i) {
        for (int8_t e = 6-offsetBit*2; e >= 0; e = e - 2)
            str[pos++] = (binaryStrAt[i] >> e) & 3;
        offsetBit = 0;
    }
    for (uint8_t i = 0; i < res % 4; ++i)
        str[pos++] = (binaryStrAt[k/4+res/4] >> (6-i*2)) & 3;
    return str;
}

inline uint8_t* bitDecodeToArr(uint8_t *binaryStr, uint8_t *str, uint64_t offset, uint64_t k) {
    
    uint8_t *binaryStrAt = &binaryStr[offset/4];
    uint8_t offsetBit = offset % 4;
    
    uint8_t res = offsetBit + k % 4;
    uint64_t pos = 0;
    for (uint64_t i = 0; i < k/4+res/4; ++i) {
        for (int8_t e = 6-offsetBit*2; e >= 0; e = e - 2)
            str[pos++] = (binaryStrAt[i] >> e) & 3;
        offsetBit = 0;
    }
    for (uint8_t i = 0; i < res % 4; ++i)
        str[pos++] = (binaryStrAt[k/4+res/4] >> (6-i*2)) & 3;
    return str;
}

inline std::string bitDecodeToStr(uint8_t *binaryStr, uint8_t offset, uint64_t k) {
    
    uint8_t *binaryStrAt = &binaryStr[offset/4];
    uint8_t offsetBit = offset % 4;
    
    std::string str;
    uint8_t res = offsetBit + k % 4;

    for (uint64_t i = 0; i < k/4+res/4; ++i) {
        for (int8_t e = 6-offsetBit*2; e >= 0; e = e - 2)
            str.push_back(itoc[(binaryStrAt[i] >> e) & 3]);
        offsetBit = 0;
    }
    for (uint8_t i = 0; i < res % 4; ++i)
        str.push_back(itoc[(binaryStrAt[k/4+res/4] >> (6-i*2)) & 3]);
    
    return str;
}

inline uint16_t bitHash(uint8_t *binaryStr, uint64_t offset, uint32_t k) {
    
    uint8_t *binaryStrAt = &binaryStr[offset/4], shift = (offset % 4)*2;
    uint16_t hash = 0;
    hash |= ((uint16_t)binaryStrAt[0] << (CHAR_BIT + shift)) + (binaryStrAt[1] << shift) + (binaryStrAt[2] >> (8-shift));
    return hash >> (8-k)*2;
}


#endif /* BIT_PACKING_H */
