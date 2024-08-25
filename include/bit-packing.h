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
struct Buf2bit : Buf<TYPE> { // 2-bit specialization of a buffer
    
    using Buf<TYPE>::Buf;
    
    Buf2bit(uint64_t size) : Buf<TYPE>((size/4+(size % 4 != 0))) { // constructor to allocate size bytes
        alloc += this->size*sizeof(TYPE);
        this->pos = size;
    }
    
    Buf2bit(std::string str) : Buf<TYPE>((str.size()/4+(str.size() % 4 != 0))) { // build 2-bit string from std::string
        
        this->pos = str.size();
        for (uint64_t i = 0; i<this->pos; ++i)
            this->seq[i / 4] |= ctoi[(unsigned char)str[i]] << (6 - (i % 4) * 2); // 2-bit packing base packing
    }
    
    uint64_t length() const { // size in characters of the 2-bit string
        return this->pos;
    }
    
    inline uint8_t at(uint64_t index) const { // return char at pos
        return (this->seq[index/4] >> (6-(index%4)*2)) & 3;
    }
    
    inline void assign(uint64_t p, char base) { // assign base ACGT
        this->seq[p/4] |= ctoi[(unsigned char)base] << (6 - (p%4) * 2);
    }
    
    inline void assign(uint64_t p, uint8_t base) { // assign base 0123
        this->seq[p/4] |= base << (6 - (p%4) * 2);
    }
    
    inline std::string substr(uint64_t offset, uint64_t k) const { // substring 2-bit string and return std::string
        
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
    
    inline std::string toString() const { // overload to use default values and convert the whole string
        return substr(0, this->pos);
    }
    
    inline void newSize(uint64_t add) { // for backward compatibility we make a new function to increase the size of the container
        
        if (this->pos + add > this->size) {
            
            uint64_t newSize = (this->pos + add)*2;
            
            alloc += newSize*sizeof(TYPE);
            TYPE* seqNew = new TYPE[newSize](); // parentheses ensure bits are set to 0
            
            memcpy(seqNew, this->seq, this->size*sizeof(TYPE));
            
            delete[] this->seq;
            freed += this->size*sizeof(TYPE);
            this->size = newSize;
            this->seq = seqNew;
        }
    }
    
    inline Buf2bit& append(const Buf2bit& buf2bit) { // append 2-bit string to 2-bit string
        
        uint64_t lenB = buf2bit.length();
        this->newSize((lenB/4+(lenB % 4 != 0))); // increase size of the container if needed
        for (uint64_t i = 0; i<lenB; ++i)
            this->assign(this->pos++,buf2bit.at(i));
        
        return *this;
    }
    
    inline uint8_t* toArray(uint64_t offset, uint64_t k) { // convert 2-bit string to uint8_t array
        
        uint8_t *str = new uint8_t[k];
        uint8_t *binaryStrAt = &this->seq[offset/4];
        uint8_t offsetBit = offset % 4;
        
        uint8_t res = offsetBit + k % 4;
        uint64_t p = 0;
        for (uint64_t i = 0; i < k/4+res/4; ++i) {
            for (int8_t e = 6-offsetBit*2; e >= 0; e = e - 2)
                str[p++] = (binaryStrAt[i] >> e) & 3;
            offsetBit = 0;
        }
        for (uint8_t i = 0; i < res % 4; ++i)
            str[p++] = (binaryStrAt[k/4+res/4] >> (6-i*2)) & 3;
        return str;
    }

    inline uint16_t bitHash(uint64_t offset, uint32_t k) { // convert position to hash
        
        uint8_t *binaryStrAt = &this->seq[offset/4], shift = (offset % 4)*2;
        uint16_t hash = 0;
        hash |= ((uint16_t)binaryStrAt[0] << (CHAR_BIT + shift)) + (binaryStrAt[1] << shift) + (binaryStrAt[2] >> (8-shift));
        return hash >> (8-k)*2;
    }
    
};

static inline Buf2bit<> revCom(Buf2bit<> &seq) { // reverse complement
    
    uint64_t len = seq.length();
    Buf2bit<> rc(len);
    
    for (uint64_t i = 0; i<len; ++i)
        rc.assign(i,static_cast<uint8_t>(3-seq.at(len-i-1)));
    return rc;
}

static inline uint16_t revCom(uint16_t hash, uint8_t k) { // reverse complement
    
    hash = hash ^ ((1 << (k*2)) - 1);
    uint16_t rc = 0;

    int8_t max = k*2-2;
    for (int8_t i = 0; i < k/2; ++i) {
        int8_t e = max-2*i;
        rc |= (hash << (e-2*i)) & (3 << e);
        rc |= (hash >> (e-2*i)) & (3 << (max-e));
    }
    if (k%2) // uneven k
        rc |= (hash & (3 << (k-1)));
    
    return rc;
}

template <typename TYPE = uint8_t>
struct Buf1bit : Buf<TYPE> { // 1-bit specialization of a buffer
    
    using Buf<TYPE>::Buf;
    
    Buf1bit(uint64_t size) : Buf<TYPE>((size/8+(size % 8 != 0))) { // constructor to allocate size bytes
        alloc += this->size*sizeof(TYPE);
        this->pos = size;
    }
    
    uint64_t length() const { // size in characters of the 2-bit string
        return this->pos;
    }
    
    inline uint8_t at(uint64_t index) const { // return char at pos
        return (this->seq[index/8] >> (7-index%8)) & 1;
    }
    
    inline void assign(uint64_t p) { // assign 1
        this->seq[p/8] |= 1 << (7 - p%8);
    }
    
    inline std::string substr(uint64_t offset, uint64_t k) { // substring 1-bit string and return std::string
        
        uint8_t *binaryStrAt = &this->seq[offset/8];
        uint8_t offsetBit = offset % 8;
        
        std::string str;
        uint8_t res = offsetBit + k % 8;
        for (uint64_t i = 0; i < k/8+res/8; ++i) {
            for (int8_t e = 7-offsetBit; e >= 0; --e)
                str.push_back('0' + ((binaryStrAt[i] >> e) & 1));
            offsetBit = 0;
        }
        for (uint8_t i = 0; i < res % 8; ++i)
            str.push_back('0' + ((binaryStrAt[k/8+res/8] >> (7-i)) & 1));
        return str;
    }
    
    inline std::string toString() { // overload to use default values and convert the whole string
        return substr(0, this->pos);
    }
    
    inline void newSize(uint64_t add) { // for backward compatibility we make a new function to increase the size of the container
        
        if (this->pos + add > this->size) {
            
            uint64_t newSize = (this->pos + add)*2;
            
            alloc += newSize*sizeof(TYPE);
            TYPE* seqNew = new TYPE[newSize](); // parentheses ensure bits are set to 0
            
            memcpy(seqNew, this->seq, this->size*sizeof(TYPE));
            
            delete[] this->seq;
            freed += this->size*sizeof(TYPE);
            this->size = newSize;
            this->seq = seqNew;
        }
    }
    
    inline Buf1bit& append(const Buf1bit& buf1bit) { // append 2-bit string to 2-bit string
        
        uint64_t lenB = buf1bit.length();
        this->newSize((lenB/8+(lenB % 8 != 0))); // increase size of the container if needed
        for (uint64_t i = 0; i<lenB; ++i) {
            if(buf1bit.at(i))
                this->assign(this->pos+i);
        }
        this->pos += lenB;
        
        return *this;
    }
};

#endif /* BIT_PACKING_H */
