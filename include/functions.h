//
//  functions.h
//
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <memory>
#include <chrono>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <string>
#include <sys/stat.h>
#include <fstream>
#include <numeric>
#include <tuple>

#ifdef EVP
#include <openssl/opensslv.h>
#include <openssl/evp.h>
#endif

#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>

#ifdef _WIN32
#include <direct.h>
#endif

#include <parallel-hashmap/phmap.h>
#include "global.h"
#include "bed.h"
#include "struct.h"

//functions

static inline std::istream& ignore(std::istream& is, char dlm) {
    
    std::ios_base::iostate err = std::ios_base::goodbit;
    std::streamsize extr = 0;
    int i = 0;
    
    while (true)
    {
        
        i = is.rdbuf()->sbumpc();
        if (i == EOF) {
            err |= std::ios_base::eofbit;
            break;
        }
        ++extr;
        if (i == dlm)
            break;
        
    }
    
    if (extr == 0)
        err |= std::ios_base::failbit;
    is.setstate(err);
    
    return is;
    
}

static inline std::istream& getKmers(std::istream& is, std::string& str, int batchSize) { // a generic function extracting concatenated reads from FASTQ
    
    str.clear();
    str.reserve(batchSize);
    std::ios_base::iostate err = std::ios_base::goodbit;
    std::streamsize extr = 0;
    int i = 0;
    
    while (batchSize > 0) {
        
        ignore(is, '\n');
        if (is.rdstate() != std::ios_base::goodbit)
            break;
        
        while (true) {
            
            i = is.rdbuf()->sbumpc();
            if (i == EOF) {
                err |= std::ios_base::eofbit;
                break;
            }
            ++extr;
            --batchSize;
            
            if (i == '\n')
                break;
            
            str.push_back(i);
            if (str.size() == str.max_size()) {
                err |= std::ios_base::failbit;
                break;
            }
            
        }

        ignore(is, '\n');
        if (is.rdstate() != std::ios_base::goodbit)
            break;
        
        ignore(is, '\n');
        if (is.rdstate() != std::ios_base::goodbit)
            break;
        
        str.push_back('N');
        
    }
    
    if (extr == 0)
        err |= std::ios_base::failbit;
    is.setstate(err);
    
    return is;
    
}

static inline std::istream& getline(std::istream& is, std::string& str) {
    
    str.clear();
    std::ios_base::iostate err = std::ios_base::goodbit;
    std::streamsize extr = 0;
    int i = 0;
    
    while (true)
    {
        i = is.rdbuf()->sbumpc();
        if (i == EOF) {
            err |= std::ios_base::eofbit;
            break;
        }
        ++extr;
        if (i == '\n')
            break;
        str.push_back(i);
        if (str.size() == str.max_size()) {
            err |= std::ios_base::failbit;
            break;
        }
    }
    
    if (extr == 0)
        err |= std::ios_base::failbit;
    is.setstate(err);
    
    return is;
    
}

static inline std::istream& getline(std::istream& is, std::string& str, char dlm) {
    
    str.clear();
    std::ios_base::iostate err = std::ios_base::goodbit;
    std::streamsize extr = 0;
    int i = 0;
    
    while (true)
    {
        
        i = is.rdbuf()->sbumpc();
        if (i == EOF) {
            err |= std::ios_base::eofbit;
            break;
        }
        ++extr;
        if (i == dlm)
            break;
        if(i == '\n')
            continue;
        str.push_back(i);
        if (str.size() == str.max_size()) {
            err |= std::ios_base::failbit;
            break;
        }
        
    }
    
    if (extr == 0)
        err |= std::ios_base::failbit;
    is.setstate(err);
    
    return is;
    
}

static inline double elapsedTime(){ // compute runtime in verbose mode
    
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
    start = std::chrono::high_resolution_clock::now();
    
    return elapsed.count();
    
}

static inline bool checkTag(const char tag1[2], std::string tag2) {
    return tag1 == tag2;
}

static inline bool isInt(const std::string &str) {
	if (str.empty()) return false;

	size_t start = 0;
	if (str[0] == '-' || str[0] == '+') {
		if (str.size() == 1) return false; // Only a sign, no digits
		start = 1;
	}

	return str.find_first_not_of("0123456789", start) == std::string::npos;
}

static inline double gfa_round(double d, uint32_t to=2) {
    if(std::isnan(d)) return NAN;
    unsigned int n=1;
    for(; to>0; to--) n *= 10;
    return std::round(d*n)/n;
}

static inline std::vector<unsigned int> intervalSizes(std::vector<unsigned int> &intervalVec){ // compute sizes of a vector of intervals to be read as paired coordinates in 0-bed format
    
    std::vector<unsigned int> intervalVecLens;
    intervalVecLens.reserve(200);
    
    if (!intervalVec.empty()){
        std::vector<unsigned int>::const_iterator end = intervalVec.cend();
        
        for (std::vector<unsigned int>::const_iterator it = intervalVec.cbegin(); it != end;) {
            
            intervalVecLens.push_back(*(it+1) - *it); // compute size of the interval
            
            it = it + 2; // jump to next pair
            
        }
        
    }
    
    return intervalVecLens;
    
}

static inline std::string output(std::string output){ // use tab delimiter if tabular flag is true
    
    if (tabular_flag) {
        output = output + "\t";
        
    }else{
        
        output = output + ": ";
        
    }
    
    return output;
}

static inline bool isDash(char * optarg) { // check if user input is dash (substitute of input from pipe)
    
    return (strcmp(optarg, "-") == 0) ? true : false;
    
}

static inline bool ifFileExists(const char * optarg) { // check if file exists
    
    if (!access (optarg, F_OK)) {
        return optarg;
    }else{
        std::cout<<"Error - file does not exist: "<<optarg<<std::endl;
        exit(1);
    }
}

static inline void textWrap(std::string input, std::ostream& output, uint32_t width) { // generic text wrapper (useful for fasta output)
    
    std::string tmp;
    char cur = '\0';
    uint32_t i = 0;
    
    std::stringstream ss(input);
    
    while (ss.get(cur)) {
        if (++i == width+1) {

            output << tmp << '\n';
            i = tmp.length();
            tmp.clear();
            
        }else{
            
            output << tmp;
            tmp.clear();
            
        }
        
        tmp += cur;
        
    }
    
    output << tmp << '\n';
    tmp.clear();
    
}

static inline std::string rmFileExt(const std::string path) { // utility to strip file extension from file
    if (path == "." || path == "..")
        return path;

    size_t pos = path.find_last_of("\\/.");
    if (pos != std::string::npos && path[pos] == '.')
        return path.substr(0, pos);

    return path;
}

static inline std::string stripKnownExt(const std::string& fname, const std::string& ext) {
	if (fname.size() >= ext.size() + 1 &&  // +1 ensures something before ext
		fname.compare(fname.size() - ext.size(), ext.size(), ext) == 0) {
		return fname.substr(0, fname.size() - ext.size() - 1);
	}
	return fname; // no change if no match
}

static inline std::string getFileExt(std::string fileName) { // utility to get file extension

    if(fileName.find_last_of(".") != std::string::npos) {
        
        if(fileName.substr(fileName.find_last_of(".")+1) == "gz") {
            fileName = rmFileExt(fileName);
            return getFileExt(fileName) + ".gz";
        }
        return fileName.substr(fileName.find_last_of(".")+1);
    }
    return "";
}

static inline std::string getFileName(std::string path) { // utility to get file extension
    return path.substr(path.find_last_of("/\\") + 1);
}

static inline std::string revCom(std::string seq) { // reverse complement
    auto lambda = [](const char c) {
        switch (c) {
        case '*':
            return '*';
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        case 'a':
            return 't';
        case 'g':
            return 'c';
        case 'c':
            return 'g';
        case 't':
            return 'a';
        case 'N':
        case 'n':
        case 'X':
        case 'x':
            return c;
        default:
            throw std::domain_error("Invalid nucleotide.");
        }
    };

    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);
    reverse(seq.begin(), seq.end());
    return seq;
}

static inline char revCom(char c) { // reverse complement
    auto lambda = [](const char c) {
        switch (c) {
        case '*':
            return '*';
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        case 'a':
            return 't';
        case 'g':
            return 'c';
        case 'c':
            return 'g';
        case 't':
            return 'a';
        case 'N':
        case 'n':
        case 'X':
        case 'x':
            return c;
        default:
            throw std::domain_error("Invalid nucleotide.");
        }
    };
    return lambda(c);
}

static inline std::string rev(std::string seq) { // reverse string
    
    reverse(seq.begin(), seq.end());
 
    return seq;
    
}

static inline std::vector<std::string> readDelimited(std::string line, std::string delimiter, std::string skipLine = "") { // read line delimited by specific character, optionally skip lines starting with specific string

    std::vector<std::string> arguments;
    size_t pos = 0;
    
    if (skipLine != "" && line.substr(0, skipLine.size()) == skipLine)
        return arguments;

    while ((pos = line.find(delimiter)) != std::string::npos) {
        arguments.push_back(line.substr(0, pos));
        line.erase(0, pos + delimiter.length());
    }
    arguments.push_back(line); // last column
    return arguments;
    
}

static inline std::vector<std::string> readDelimitedArr(std::string line, std::vector<char> delimiters, std::string skipLine = "", bool keepDelimiter = false) { // read line delimited by specific character, optionally skip lines starting with specific string
    
    auto is_delimiter = [delimiters](char c){
        
        bool cond = false;
        
        for (char delimiter : delimiters) {
            
            cond = (c == delimiter ? true : false);
            
            if (cond)
                break;
            
        }
        
        return cond;
        
    };

    std::vector<std::string> arguments;

    size_t pos = 0;
    
    if (skipLine != "" && line.substr(0, skipLine.size()) == skipLine) {
        
        return arguments;
        
    }
    
    auto it = begin(line);

    while ((it = std::find_if(begin(line), end(line), is_delimiter)) != end(line)) {
        
        pos = std::distance(begin(line), it);
        
        arguments.push_back(line.substr(0, pos + (keepDelimiter ? 1 : 0)));
        
        line.erase(0, pos + (keepDelimiter ? 1 : 0));
            
    }
    
    arguments.push_back(line); // last column
        
    return arguments;
    
}

static inline bool isNumber(const std::string& str)
{
    for (char const &c : str) {
        if (std::isdigit(c) == 0) return false;
    }
    return true;
}

static inline void revComPathComponents(std::vector<PathComponent>& pathComponents) {
    
    for (PathComponent& component : pathComponents) {
        
        if (component.orientation != '0') {
        
            component.orientation = (component.orientation == '+' ? '-' : '+');
        
        }
        
    }
    
    std::reverse(pathComponents.begin(), pathComponents.end());
    
}

// bed coords are bed coords of compressed sequence
static inline void homopolymerCompress(std::string *sequence, std::vector<std::pair<uint64_t, uint64_t>> &bedCoords, unsigned int cutoff) {
    unsigned int index=0, length, new_length=0;

    auto lambda = [&length, &index, &bedCoords, &sequence, &new_length, &cutoff](int i){
        length = i-index;
        if(length > cutoff) {
            bedCoords.push_back({new_length, new_length+length});
        }
        int num = length > cutoff ? 1 : length;
        memset(&((*sequence)[new_length]), (*sequence)[index], num);
        new_length += num;
    };

    for(unsigned int i=1; i<sequence->length(); ++i) {
        if((*sequence)[i] == (*sequence)[index]) continue;
        lambda(i);
        index = i;
    }
    lambda(sequence->length());
    sequence->resize(new_length);
}

// bed coords are bed coords of compressed sequence
static inline void homopolymerDecompress(std::string *sequence, const std::vector<std::pair<uint64_t, uint64_t>> &bedCoords) {
    std::string ret="";
    ret.reserve(sequence->length()*2); // random guess for final sequence length to reduce resizes
    for(unsigned int i=0, ci=0, len; i<sequence->length(); ++i) {
        if(ci < bedCoords.size() && i == bedCoords[ci].first) {
            len = bedCoords[ci].second - bedCoords[ci].first;
            ++ci;
        } else {
            len = 1;
        }
        ret += std::string(len, (*sequence)[i]);
    }
    ret.shrink_to_fit();
    *sequence = ret;
}

static inline unsigned int homopolymerRunsCount(const std::string &sequence, unsigned int threshhold) {
    unsigned int runs = 0;
    unsigned int currentRun=0;
    char prev=0;
    for(const char &curr : sequence) {
        if(prev != curr) {
            if(currentRun >= threshhold) {
                ++runs;
            }
            currentRun = 0;
        }
        else {
            ++currentRun;
        }
        prev = curr;
    }
    if(currentRun >= threshhold) { // loop wont catch run at end of sequence
        ++runs;
    }
    return runs;
}

// bed coords of uncompressed sequence
static inline void homopolymerBedCoords(std::string *sequence, std::vector<std::pair<unsigned int, unsigned int>> &bedCoords, unsigned int cutoff) {
    int index = 0;
    for(unsigned int i=1; i < sequence->size(); ++i) {
        if((*sequence)[i] == (*sequence)[i-1]) continue;
        if(i-index > cutoff) {
            bedCoords.push_back({index, i});
        }
        index = i;
    }
    if((*sequence)[sequence->size()-1] == (*sequence)[sequence->size()-2] && sequence->size()-index > cutoff) {
        bedCoords.push_back({index, sequence->size()});
    }
}

static inline void computeNstars(std::vector<uint64_t>& lens, // compute N/L* statistics, vector of all lengths, returns sorted vector of lengths
                   std::vector<uint64_t>& Nstars,      std::vector<unsigned int>& Lstars, // required arguments are passed by reference
                   std::vector<uint64_t>* NGstars = NULL, std::vector<unsigned int>* LGstars = NULL, uint64_t gSize = 0) { // optional arguments are passed by pointer
    
    sort(lens.begin(), lens.end(), std::greater<uint64_t>()); // sort lengths Z-A
    uint64_t sum = 0, totLen = 0;
    
    for(std::vector<uint64_t>::iterator it = lens.begin(); it != lens.end(); ++it) // find total length
        totLen += *it;
    
    short int N = 1, NG = 1;
    
    for(unsigned int i = 0; i < lens.size(); i++) { // for each length
        
        sum += lens[i]; // increase sum
        while (sum >= ((double) totLen / 10 * N) && N<= 10) { // conditionally add length.at or pos to each N/L* bin
            
            Nstars[N-1] = lens[i];
            Lstars[N-1] = i + 1;
            N = N + 1;
        }
        while (gSize > 0 && (sum >= ((double) gSize / 10 * NG)) && NG<= 10) { // if not computing gap statistics repeat also for NG/LG* statistics
            
            (*NGstars)[NG-1] = lens[i];
            (*LGstars)[NG-1] = i + 1;
            NG = NG + 1;
        }
    }
}

static inline void rmChrFromStr(std::string &str, const char* charsToRemove) {
   for (unsigned int i = 0; i < strlen(charsToRemove); ++i ) {
      str.erase(std::remove(str.begin(), str.end(), charsToRemove[i]), str.end());
   }
}

static inline void make_dir(const char* name) {
#ifdef _WIN32
    _mkdir(name);
#else
    mkdir(name, 0777);
#endif
}

static inline void rm_dir(const char* name) {
#ifdef _WIN32
    _rmdir(name);
#else
    rmdir(name);
#endif
}

static inline bool fileExists(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

static inline unsigned int fileCount(const char *dir) {
    struct dirent *dp;
    DIR *fd;
    unsigned int i = 0;

    if ((fd = opendir(dir)) == NULL) {
        fprintf(stderr, "Error: can't open %s\n", dir);
        exit(1);
    }
    while ((dp = readdir(fd)) != NULL) {
        if (!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, ".."))
            continue;    /* skip self and parent */
        ++i;
    }
    closedir(fd);
    return i;
}

static inline uint64_t fileSize(std::string path) {
    
    ifFileExists(path.c_str());
    std::ifstream file(path, std::ios::binary);
    const auto begin = file.tellg();
    file.seekg(0, std::ios::end);
    const auto end = file.tellg();
    return end-begin;
    
}

static inline std::vector<uint32_t> sortedIndex(std::vector<uint64_t> vec, bool largest) {
    
    std::vector<uint32_t> idx(vec.size());
    std::iota(idx.begin(),idx.end(),0);
    if (largest)
        stable_sort(idx.begin(), idx.end(), [&vec](size_t i,size_t j){return vec[i]>vec[j];} );
    else
        stable_sort(idx.begin(), idx.end(), [&vec](size_t i,size_t j){return vec[i]<vec[j];} );
    
    return idx;
    
}

static inline void unmaskSequence(std::string &sequence) {
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
}

static inline void eraseChar(std::string& input, char rmChar) {
    input.erase(std::remove(input.begin(), input.end(), rmChar), input.end());
}

static inline std::tuple<std::string, uint64_t, uint64_t> parseCoordinate(std::string input) {
    
    std::string header, cBegin, cEnd; // the header for coordinates provided as positional argument
    uint64_t cBeginNumeric = 0, cEndNumeric = 0;
    
    reverse(input.begin(), input.end()); // we work our way from the end
    
    cBegin = input.substr(input.find('-') + 1, input.find(':') - input.find('-') - 1);
    cEnd = input.substr(0, input.find('-'));
    
    if(isNumber(cEnd) && isNumber(cBegin)) { // prevent headers with : - characters to break the extraction
        
        header = input.substr(input.find(':') + 1, input.size());
        reverse(header.begin(), header.end());
        
        reverse(cBegin.begin(), cBegin.end());
        reverse(cEnd.begin(), cEnd.end());
        
        cBeginNumeric = stoull(cBegin);
        cEndNumeric = stoull(cEnd);
        
    }else{
        header = input;
        reverse(header.begin(), header.end());
    }
    
    return std::make_tuple(header, cBeginNumeric, cEndNumeric);
}

static inline uint64_t parseCigar(std::string cigar) { // only works with M (identity)
	uint64_t pos = cigar.find_first_of('M');
	return stoi(cigar.substr(0, pos));
}

#ifdef EVP
static inline bool computeMd5(const std::string file, std::string &md5) {
    unsigned char md_value[EVP_MAX_MD_SIZE];
    unsigned int  md_len;

    EVP_MD_CTX*   context = EVP_MD_CTX_new();
    const EVP_MD* md = EVP_md5();

#if OPENSSL_VERSION_NUMBER >= 0x30000000L
    EVP_DigestInit_ex2(context, md, NULL);
#else
    EVP_DigestInit_ex(context, md, NULL);
#endif

    const int bufSize = 1024;
    char buffer[bufSize];

    std::ifstream fin(file);

    while(!fin.eof()) {
        fin.read(buffer, bufSize);
        std::streamsize s=fin.gcount();
        EVP_DigestUpdate(context, buffer, s);
    }

    EVP_DigestFinal_ex(context, md_value, &md_len);
    EVP_MD_CTX_free(context);

    char *md_value_buf = new char[md_len*2+1]; // two characters per digit + null termination

    for (unsigned int i = 0 ; i < md_len ; ++i)
        snprintf(md_value_buf+i*2, 3, "%02x", md_value[i]);

    md5 = std::string(md_value_buf);
    delete[] md_value_buf;
    return true;
}
#endif

#endif /* FUNCTIONS_H */
