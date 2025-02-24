#ifndef SAK_H
#define SAK_H

#include <stack>

struct Instruction {
    std::string action,
                path1,
                path2,
                path3,
                contig1,
                contig2,
                contig3,
                gap1,
                comment1,
                comment2,
                gHeader,
                edge1;
    
    int compressThreshhold;
    char pId1Or, pId2Or, sId1Or, sId2Or;
    uint64_t gUId = 0, dist = 0, start1 = 0, end1 = 0, start2 = 0, end2 = 0;
};

class SAK;

struct Funcs {
public:
    void (SAK::*read)(Instruction&, std::vector<std::string>&);
    bool (SAK::*exec)(InSequences&, Instruction&);
};

class SAK { // the swiss army knife
public:

    void readJoin(Instruction &instruction, std::vector<std::string> &arguments) {
        size_t pos = 0, pos2 = 0, pos3 = 0;
        
        instruction.pId1Or = arguments[1].back(); // get sequence orientation in the gap
        arguments[1].pop_back(); // remove sequence orientation in the gap
        
        if(arguments[1].back() == ')') {
            pos = arguments[1].find('(');
            pos2 = arguments[1].find(':');
            pos3 = arguments[1].find(')');
            
            instruction.path1 = arguments[1].substr(0, pos);
            instruction.start1 = stoi(arguments[1].substr(pos + 1, pos2 - pos - 1));
            instruction.end1 = stoi(arguments[1].substr(pos2 + 1, pos3 - pos2 - 1));
        }else{
            instruction.path1 = arguments[1];
        }
    
        instruction.pId2Or = arguments[2].back(); // get sequence orientation in the gap
        arguments[2].pop_back(); // remove sequence orientation in the gap
        
        if(arguments[2].back() == ')') {
            pos = arguments[2].find('(');
            pos2 = arguments[2].find(':');
            pos3 = arguments[2].find(')');
            
            instruction.path2 = arguments[2].substr(0, pos);
            instruction.start2 = stoi(arguments[2].substr(pos + 1, pos2 - pos - 1));
            instruction.end2 = stoi(arguments[2].substr(pos2 + 1, pos3 - pos2 - 1));
        }else{
            instruction.path2 = arguments[2];
        }
        
        instruction.dist = stoi(arguments[3]);
        instruction.gHeader = arguments[4];
        instruction.path3 = arguments[5];
    }
    bool join(InSequences& inSequences, Instruction &instruction) { // joins two sequences via a gap based on instruction
        inSequences.joinPaths(instruction.path3, inSequences.headersToIds[instruction.path1], inSequences.headersToIds[instruction.path2], instruction.gHeader, instruction.gUId, instruction.pId1Or, instruction.pId2Or, instruction.dist, instruction.start1, instruction.end1, instruction.start2, instruction.end2); // generate a new path by joining the paths that contain the two segments
        return true;
    }
    void readSplit(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
        instruction.contig2 = arguments[2];
        instruction.path1 = arguments[3];
        instruction.path2 = arguments[4];
    }
    bool split(InSequences& inSequences, Instruction &instruction) { // splits two sequences removing the gap in between based on instruction
        std::vector<unsigned int> guIds = inSequences.removeGaps(&instruction.contig1, &instruction.contig2); // remove the gap
        inSequences.splitPath(guIds[0], instruction.path1, instruction.path2); // generate two new the paths splitting the original path
        return true;
    }
    void readExcise(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];

        if (arguments[2] != "") {
            instruction.dist = stoi(arguments[2]);
        }else{
            instruction.dist = 0;
        }
        instruction.gHeader = arguments[3];
    }
    bool excise(InSequences& inSequences, Instruction &instruction) { // excises a sequence, removing also edges if present and optionally adding a gap
        std::vector<InGap> oldGaps = inSequences.getGap(&instruction.contig1); // get the gaps associated with contig1
        InGap gap;

        gap.newGap(inSequences.uId.get(), oldGaps[0].getsId1(), oldGaps[1].getsId2(), oldGaps[0].getsId1Or(), oldGaps[1].getsId2Or(), instruction.dist, instruction.gHeader); // define the new gap
        
        inSequences.insertHash(instruction.gHeader, inSequences.uId.get());
        
        inSequences.uId.next();
        
        inSequences.addGap(gap); // introduce the new gap
        inSequences.removeGaps(&instruction.contig1);
        inSequences.removeSegmentInPath(inSequences.headersToIds[instruction.contig1], gap); // removes the segment from the path
        return true;
    }
    void readRemove(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
    }
    bool remove(InSequences& inSequences, Instruction &instruction) { // removes a segment
        inSequences.removePathsFromSegment(inSequences.headersToIds[instruction.contig1]); // remove the paths involving contig1
        inSequences.flagDeletedSegment(&instruction.contig1); // flag segment as deleted
        return true;
    }
    void readExclude(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.path1 = arguments[1];
    }
    bool exclude(InSequences& inSequences, Instruction &instruction) { // removes a segment
        inSequences.removePath(inSequences.headersToIds[instruction.path1], true); // remove the path and its components
        return true;
    }
    void readErase(Instruction &instruction, std::vector<std::string> &arguments) {  
        size_t pos1 = 0, pos2 = 0;
        pos1 = arguments[1].find(":");
        pos2 = arguments[1].find("-");
        instruction.contig1 = arguments[1].substr(0, pos1);
        instruction.start1 = stoi(arguments[1].substr(pos1+1, pos2));
        instruction.end1 = stoi(arguments[1].substr(pos2+1, arguments[1].size()+1));
    }
    bool erase(InSequences& inSequences, Instruction &instruction) { // erases a portion of sequence
        
        lg.verbose("Erasing from " + instruction.contig1 + " coordinates " + std::to_string(instruction.start1) + ":" + std::to_string(instruction.end1));
        
        unsigned int uId = inSequences.headersToIds[instruction.contig1], sIdx = 0;
        auto sId = find_if(inSequences.inSegments.begin(), inSequences.inSegments.end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node uId, find it

        if (sId != inSequences.inSegments.end()) {sIdx = std::distance(inSequences.inSegments.begin(), sId);} // gives us the segment index
        
        inSequences.inSegments[sIdx]->trimSegment(instruction.start1, instruction.end1); // trim segment
        inSequences.inSegments[sIdx]->updateSegmentCounts(instruction.start1, instruction.end1);
        return true;
    }
    void readRvcp(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
    }
    bool rvcp(InSequences& inSequences, Instruction &instruction) { // reverse complement sequence
        
        unsigned int uId = inSequences.headersToIds[instruction.contig1], sIdx = 0;
        
        auto sId = find_if(inSequences.inSegments.begin(), inSequences.inSegments.end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node uId, find it
        
        auto pId = find_if(inSequences.inPaths.begin(), inSequences.inPaths.end(), [uId](InPath obj) {return obj.getpUId() == uId;}); // given a node uId, find it

        if (sId != inSequences.inSegments.end()) {
            
            sIdx = std::distance(inSequences.inSegments.begin(), sId); // gives us the segment index
            
            inSequences.inSegments[sIdx]->rvcpSegment(); // rvcp segment
            
        }else if (pId != inSequences.inPaths.end()){
            
            inSequences.revComPath(uId);
            
        }else{
            
            fprintf(stderr, "Error: couldn't find (%s). Terminating.\n", instruction.contig1.c_str());
            exit(1);
            
        }
        
        return true;
    }
    void readInvert(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
    }
    bool invert(InSequences& inSequences, Instruction &instruction) { // invert sequence
        unsigned int uId = inSequences.headersToIds[instruction.contig1], sIdx = 0;
        auto sId = find_if(inSequences.inSegments.begin(), inSequences.inSegments.end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node uId, find it

        if (sId != inSequences.inSegments.end()) {sIdx = std::distance(inSequences.inSegments.begin(), sId);} // gives us the segment index
        inSequences.inSegments[sIdx]->invertSegment(); // invert segment
        return true;
    }
    void readCompress(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
        instruction.compressThreshhold = stoi(arguments[2]);
    }
    std::string *headerToSequence(const std::string &header, InSequences &inSequences, const Instruction &instruction) {
        auto got = inSequences.headersToIds.find(header); // get the headers to uIds table to look for the header
        if (got == inSequences.headersToIds.end()) { // this is the first time we see this path name
            fprintf(stderr, "Error: couldn't find (%s). Terminating.\n", instruction.contig1.c_str());
            exit(1);
        }
        int id = got->second;
        return (inSequences.getInSegment(id)->inSequence);
    }
    bool compress(InSequences &inSequences, Instruction &instruction) {
        std::string *sequence = headerToSequence(instruction.contig1, inSequences, instruction);
        std::vector<std::pair<uint64_t, uint64_t>> bedCoords;
        homopolymerCompress(sequence, bedCoords, instruction.compressThreshhold);
        compressStack.push(bedCoords);
        return true;
    }
    void readDecompress(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
    }
    bool decompress(InSequences &inSequences, Instruction &instruction) {
        std::string *sequence = headerToSequence(instruction.contig1, inSequences, instruction);
        std::vector<std::pair<uint64_t, uint64_t>> bedCoords = compressStack.top();
        compressStack.pop();
        homopolymerDecompress(sequence, bedCoords);
            
        return true;
    }
    void readResize(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.gap1 = arguments[1];
        instruction.dist = stoi(arguments[2]);
    }
    bool resize(InSequences &inSequences, Instruction &instruction) {
        inSequences.resizeGap(instruction.gap1, instruction.dist);
        
        return true;
    }
    void readMask(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.path1 = arguments[1];
        instruction.start1 = stoi(arguments[2]);
        instruction.end1 = stoi(arguments[3]);
        if (arguments[4] != "") {
            instruction.dist = stoi(arguments[4]);
        }else{
            instruction.dist = instruction.end1 - instruction.start1;
        }
    }
    bool mask(InSequences &inSequences, Instruction &instruction) {
        inSequences.maskPath(instruction.path1, instruction.start1, instruction.end1, instruction.dist);
        
        return true;
    }
    void readCleave(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
        instruction.start1 = stoi(arguments[2]);
        instruction.contig2 = arguments[3];
        instruction.contig3 = arguments[4];
        if (arguments[5] != "")
            instruction.edge1 = arguments[5];
    }
    bool cleave(InSequences &inSequences, Instruction &instruction) {
        inSequences.cleaveSegment(inSequences.headersToIds[instruction.contig1], instruction.start1, instruction.contig2, instruction.contig3, instruction.edge1);
        return true;
    }
    void readRename(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.path1 = arguments[1];
        instruction.path2 = arguments[2];
    }
    bool rename(InSequences &inSequences, Instruction &instruction) {
        inSequences.renamePath(instruction.path1, instruction.path2);
        return true;
    }
    void readComment(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.path1 = arguments[1];
        instruction.comment1 = arguments[2];
    }
    bool comment(InSequences &inSequences, Instruction &instruction) {
        inSequences.updateComment(instruction.path1, instruction.comment1);
        return true;
    }
private:

    std::string sId1Header, sId2Header;
    std::stack<std::vector<std::pair<uint64_t, uint64_t>>> compressStack; // decompress uses indices of corresponding compresses, last compress = first decompress

    // unordered map from sak instruction to corresponding read and execute functions
    const phmap::flat_hash_map<std::string, Funcs> string_to_funcs {
        {"JOIN",        {&SAK::readJoin,        &SAK::join      } },
        {"SPLIT",       {&SAK::readSplit,       &SAK::split     } },
        {"EXCISE",      {&SAK::readExcise,      &SAK::excise    } },
        {"REMOVE",      {&SAK::readRemove,      &SAK::remove    } },
        {"ERASE",       {&SAK::readErase,       &SAK::erase     } },
        {"RVCP",        {&SAK::readRvcp,        &SAK::rvcp      } },
        {"INVERT",      {&SAK::readInvert,      &SAK::invert    } },
        {"COMPRESS",    {&SAK::readCompress,    &SAK::compress  } },
        {"DECOMPRESS",  {&SAK::readDecompress,  &SAK::decompress} },
        {"EXCLUDE",     {&SAK::readExclude,     &SAK::exclude   } },
        {"RESIZE",      {&SAK::readResize,      &SAK::resize    } },
        {"MASK",        {&SAK::readMask,        &SAK::mask      } },
        {"CLEAVE",      {&SAK::readCleave,      &SAK::cleave    } },
        {"RENAME",      {&SAK::readRename,      &SAK::rename    } },
        {"COMMENT",     {&SAK::readComment,     &SAK::comment   } }
    };

public:
    
    Instruction readInstruction(std::string line) {
        
        std::string delimiter = "\t";
        std::vector<std::string> arguments;
        
        size_t pos = 0;
        
        while ((pos = line.find(delimiter)) != std::string::npos) {
            
            arguments.push_back(line.substr(0, pos));
            
            line.erase(0, pos + delimiter.length());
        
        }
        
        arguments.push_back(line); // last column
        
        Instruction instruction;
        
        for (auto & c: arguments[0]) instruction.action += (char) toupper(c);
        
        if(!string_to_funcs.count(instruction.action)) {
            fprintf(stderr, "unrecognized action %s\n", instruction.action.c_str());
            exit(EXIT_FAILURE);
        }

        (this->*(string_to_funcs.at(instruction.action).read))(instruction, arguments);
        
        lg.verbose("Instruction read");
        
        return instruction;
        
    }
    
    bool executeInstruction(InSequences& inSequences, Instruction instruction) {
        if(!string_to_funcs.count(instruction.action)) {
            fprintf(stderr, "unrecognized action %s\n", instruction.action.c_str());
            exit(EXIT_FAILURE);
        }
        return (this->*(string_to_funcs.at(instruction.action).exec))(inSequences, instruction);
    }

};

#endif /* SAK_H */
