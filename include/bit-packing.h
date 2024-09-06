#ifndef BIT_PACKING_H
#define BIT_PACKING_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <bitset>
#include <set>

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

// function to convert decimal to hexadecimal
static inline std::string decToHexa(int n) {
    // ans string to store hexadecimal number
    std::string ans = "";
  
    while (n != 0) {
        // remainder variable to store remainder
        int rem = 0;
        
        // ch variable to store each character
        char ch;
        // storing remainder in rem variable.
        rem = n % 16;

        // check if temp < 10
        if (rem < 10) {
            ch = rem + 48;
        }
        else {
            ch = rem + 55;
        }
        
        // updating the ans string with the character variable
        ans += ch;
        n = n / 16;
    }
    
    // reversing the ans string to get the final result
    int i = 0, j = ans.size() - 1;
    while(i <= j)
    {
      std::swap(ans[i], ans[j]);
      i++;
      j--;
    }
    return ans;
}

template <typename TYPE = uint8_t>
struct Buf2bit : Buf<TYPE> { // 2-bit specialization of a buffer
    
    using Buf<TYPE>::Buf;
    
    Buf2bit(uint64_t size) : Buf<TYPE>((size/4+(size % 4 != 0))) { // constructor to allocate size bytes
        alloc += this->size*sizeof(TYPE);
        this->pos = size;
    }
    
    Buf2bit(std::string str) : Buf<TYPE>((str.size()/4+(str.size() % 4 != 0))) { // build 2-bit string from std::string
        
        this->pos = str.size();
        for (uint64_t i = 0; i<this->pos; ++i) {
            uint8_t base = ctoi[(unsigned char)str[i]];
            if (base < 4)
                this->seq[i / 4] |= base << (6 - (i % 4) * 2); // 2-bit packing base packing
        }
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
        for (uint8_t i = offsetBit; i < res % 4; ++i)
            str.push_back(itoc[(binaryStrAt[k/4+res/4] >> (6-i*2)) & 3]);
        
        return str;
    }
    
    inline Buf2bit substr2bit(uint64_t offset, uint64_t k) const { // substring 2-bit string and return 2-bit buffer
        
        uint8_t *binaryStrAt = &this->seq[offset/4];
        uint8_t offsetBit = offset % 4;
        
        Buf2bit str(k);
        uint8_t res = offsetBit + k % 4;
        uint64_t pos = 0;

        for (uint64_t i = 0; i < k/4+res/4; ++i) {
            for (int8_t e = 6-offsetBit*2; e >= 0; e = e - 2)
                str.assign(pos++, static_cast<uint8_t>((binaryStrAt[i] >> e) & 3));
            offsetBit = 0;
        }
        for (uint8_t i = offsetBit; i < res % 4; ++i)
            str.assign(pos++, static_cast<uint8_t>(binaryStrAt[k/4+res/4] >> (6-i*2) & 3));
        
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

    template <uint16_t (*filterFN)(uint16_t, uint8_t)>
    inline uint16_t bitHash(uint64_t offset, uint8_t k) { // convert position to hash
        
        uint8_t *binaryStrAt = &this->seq[offset/4], shift = (offset % 4)*2;
        uint16_t hash = 0;
        hash |= ((uint16_t)binaryStrAt[0] << (CHAR_BIT + shift)) + (binaryStrAt[1] << shift) + (binaryStrAt[2] >> (8-shift));
        
        hash = filterFN(hash, k);
        
        return hash >> (8-k)*2;
    }
};

// non-canonical minimizers
static inline uint16_t hashNoFilter(uint16_t hash, uint8_t k = 7) {
    (void)k;
    return hash;
}

static inline uint16_t hashNC(uint16_t hash, uint8_t k) {
    if (((hash & 0xFC00) >> 10) == 0) // AAA
        return 16385;
    
    if (((hash & 0xFC00) >> 10) == 4) // ACA
        return 16385;
    
    for (uint16_t i = 1; i<k-1; ++i) { // NAAN
        if ((hash & (15 << 2*(6-i))) == 0)
            return 16385;
    }
    return hash;
}

static inline Buf2bit<> revCom(Buf2bit<> &seq) { // reverse complement
    
    uint64_t len = seq.length();
    Buf2bit<> rc(len);
    
    for (uint64_t i = 0; i<len; ++i)
        rc.assign(i,static_cast<uint8_t>(3-seq.at(len-i-1)));
    return rc;
}

template <uint16_t (*filterFN)(uint16_t, uint8_t)>
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
    
    rc = filterFN(rc << (8-k)*2, k);
    
    return rc >> (8-k)*2;
}

template <typename TYPE = uint8_t>
struct Buf1bit : Buf<TYPE> { // 1-bit specialization of a buffer
    
    using Buf<TYPE>::Buf;
    
    Buf1bit(uint64_t size) : Buf<TYPE>((size/8+(size % 8 != 0))) { // constructor to allocate size bytes
        alloc += this->size*sizeof(TYPE);
        this->pos = size;
    }
    
    Buf1bit(std::string str) : Buf<TYPE>((str.size()/8+(str.size() % 8 != 0))) { // build 1-bit string from std::string
        
        this->pos = str.size();
        for (uint64_t i = 0; i<this->pos; ++i)
            this->seq[i / 8] |= (str[i] - '0') << (7 - i % 8); // 1-bit packing
    }
    
    uint64_t length() const { // size in characters of the 2-bit string
        return this->pos;
    }
    
    inline uint8_t at(uint64_t index) const { // return char at pos
        return (this->seq[index/8] >> (7 - index % 8)) & 1;
    }
    
    inline void assign(uint64_t p) { // assign 1
        this->seq[p/8] |= 1 << (7 - p % 8);
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
        for (uint8_t i = offsetBit; i < res % 8; ++i)
            str.push_back('0' + ((binaryStrAt[k/8+res/8] >> (7-i)) & 1));
        return str;
    }

    inline Buf1bit substr1bit(uint64_t offset, uint64_t k) { // substring 1-bit string and return buffer
        
        uint8_t *binaryStrAt = &this->seq[offset/8];
        uint8_t offsetBit = offset % 8;
        
        Buf1bit str(k);
        uint8_t res = offsetBit + k % 8;
        uint64_t pos = 0;

        for (uint64_t i = 0; i < k/8+res/8; ++i) {
            for (int8_t e = 7-offsetBit; e >= 0; --e) {
                if((binaryStrAt[i] >> e) & 1)
                    str.assign(pos);
                ++pos;
            }
            offsetBit = 0;
        }
        for (uint8_t i = offsetBit; i < res % 8; ++i) {
            if((binaryStrAt[k/8+res/8] >> (7-i)) & 1)
                str.assign(pos);
            ++pos;
        }
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
    
    inline Buf1bit& append(const Buf1bit& buf1bit) { // append 1-bit string to 1-bit string
        
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

class String2bit {
  
    Buf2bit<> array;
    Buf1bit<> mask;

public:
    String2bit(std::string str, bool generateMask = true) : array(str), mask(str.size()) { // build 2-bit string from std::string
        
        if (generateMask) {
            uint64_t strLen = str.size();
            
            for (uint64_t i = 0; i<strLen; ++i) {
                if (ctoi[(unsigned char)str[i]] > 3)
                    mask.assign(i);
            }
        }
    }
    
    String2bit(Buf2bit<> array, Buf1bit<> mask) : array(std::move(array)), mask(std::move(mask)) {}
    
    uint16_t bitHash(uint64_t offset, uint8_t k) {
        return this->array.bitHash<hashNoFilter>(offset, k);
    }
    
    std::string toString() {
        
        std::string str = array.toString();
        uint64_t strLen = str.size();
        
        for (uint64_t i = 0; i<strLen; ++i) {
            if (mask.at(i))
                str[i] = 'N';
        }
        return str;
    }
    
    String2bit substr(uint64_t offset, uint64_t k) {
        return String2bit(array.substr2bit(offset, k), mask.substr1bit(offset, k));
    }
    
    Buf2bit<>& data() {
        return this->array;
    }

    bool maskAt(uint64_t pos) {
        return mask.at(pos);
    }
    
    std::string maskToString() {
        return mask.toString();
    }
    
    uint64_t size() const { // size in characters of the 2-bit string
        return mask.pos;
    }
    
    template <uint16_t (*filterFN)(uint16_t, uint8_t)>
    Buf1bit<> minimizersToMask(uint32_t k, uint8_t s) {
        
        std::multiset<uint16_t> smers;
        uint64_t kcount = this->size()-k+1;
        Buf1bit<> minimizerMask(this->size());
        uint32_t e = 0;
        uint32_t max = std::pow(4,s)+1, minimizer = max;
        
        for (uint64_t p = 0; p<kcount; ++p) {
            
            for (uint32_t c = e; c<k; ++c) { // generate k bases if e=0 or the next if e=k-1

                if (!this->mask.at(p+c)) { // filter non ACGTacgt bases

                    if (c+1 >= s) {
                        uint16_t hash = this->array.bitHash<filterFN>(p+c+1-s, s), hashRc = revCom<filterFN>(hash, s);
                        smers.insert(hash < hashRc ? hash : hashRc);
                    }
                }
                else { // if non-canonical/N base is found
                    p = p+c; // move position
                    e = 0; // reset base counter
                    smers.clear();
                    minimizer = max;
                    break;
                }
                e = k-1; // after the first kmer we only read one char at a time
            }
            if (e == 0) // not enough bases for a kmer
                continue;
            
            if (smers.find(minimizer) == smers.end())
                minimizer = max;
            
            if (*smers.begin() < minimizer) {
                minimizerMask.assign(p);
                minimizer = *smers.begin();
            }
            uint16_t hash = this->array.bitHash<filterFN>(p, s), hashRc = revCom<filterFN>(hash, s);
            smers.erase(smers.find(hash < hashRc ? hash : hashRc));
        }
        return minimizerMask;
    }
    
    template <uint16_t (*filterFN)(uint16_t, uint8_t)>
    std::string minimizersHexString(uint32_t k, uint8_t s) {
        
        std::multiset<uint16_t> smers;
        uint64_t kcount = this->size()-k+1;
        std::string minimizerMask;
        uint32_t e = 0;
        
        for (uint64_t p = 0; p<kcount; ++p) {
            
            for (uint32_t c = e; c<k; ++c) { // generate k bases if e=0 or the next if e=k-1

                if (!this->mask.at(p+c)) { // filter non ACGTacgt bases

                    if (c+1 >= s) {
                        uint16_t hash = this->array.bitHash<filterFN>(p+c+1-s, s), hashRc = revCom<filterFN>(hash, s);
                        smers.insert(hash < hashRc ? hash : hashRc);
                    }
                }
                else { // if non-canonical/N base is found
                    p = p+c; // move position
                    e = 0; // reset base counter
                    smers.clear();
                    minimizerMask += std::string(k,'N');
                    break;
                }
                e = k-1; // after the first kmer we only read one char at a time
            }
            if (e == 0) // not enough bases for a kmer
                continue;
            
            minimizerMask += decToHexa(*smers.begin());
            uint16_t hash = this->array.bitHash<filterFN>(p, s), hashRc = revCom<filterFN>(hash, s);
            smers.erase(smers.find((hash < hashRc ? hash : hashRc)));
        }
        return minimizerMask;
    }
    
    template <uint16_t (*filterFN)(uint16_t, uint8_t)>
    uint32_t getMinimizer(uint8_t s) {
        uint64_t kcount = this->size()-s+1;
        uint32_t minimizer = std::pow(4,s)+1;
        
        for (uint64_t p = 0; p<kcount; ++p) {
            uint16_t hash = this->array.bitHash<filterFN>(p, s), hashRc = revCom<filterFN>(hash, s);
            uint16_t smer = (hash < hashRc ? hash : hashRc);
            if (smer < minimizer)
                minimizer = smer;
        }
        return minimizer;
    }
};

template <uint16_t (*filterFN)(uint16_t, uint8_t)>
class MinimizerStream {
    
    uint64_t substringStart = 0, i = 0;
    uint32_t k, gap = 0;
    String2bit &str;
    bool state = true;
    Buf1bit<> minimizerMask;
    std::multiset<uint16_t> smers;
    
public:
    
    MinimizerStream(String2bit &str, uint32_t k, uint8_t s) : k(k), str(str), minimizerMask(str.minimizersToMask<filterFN>(k,s)) {
        if (minimizerMask.length() < k)
            state = false;
    }
    
    explicit operator bool() const noexcept {
        return (state);
    }
    
    uint32_t gapSize() {
        return gap;
    }
    
    uint64_t pos() {
        return i;
    }
    
    String2bit next() {
        
        if (!state)
            return str.substr(0,0);
        
        uint64_t pos = 0, size = 0, kcount = minimizerMask.length()-k+1;
        
        for (; i < kcount; ++i) {
            
            if (minimizerMask.at(i+1)) {
                size = i-substringStart+k;
                pos = substringStart;
                substringStart = i+1;
                gap = 0;
                break;
            } else if (str.maskAt(i+1)) {
                size = i-substringStart+1;
                pos = substringStart;
                uint32_t jump = 1;
                while (str.maskAt(i+jump))
                    ++jump;
                substringStart = i+jump;
                i += jump-1;
                gap = jump-1;
                break;
            } else if (i+1 == kcount) {
                uint32_t offset = 0;
                for (uint32_t e = 0; e < k; ++e) {
                    if (str.maskAt(i+offset))
                        break;
                    ++offset;
                }
                size = i-substringStart+offset;
                pos = substringStart;
                gap = 0;
                state = false;
                break;
            }
        }
        ++i;
        return str.substr(pos, size);
    }
};
#endif /* BIT_PACKING_H */
