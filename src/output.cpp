#include <stdlib.h>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <deque>

#include <parallel-hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "log.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"

#include "output.h"

OutputStream::OutputStream(std::string file) : file(file), ofs(file), zfout(ofs), zout(std::cout) {
    
    if (file == "")
        return;
    
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal points
    
    // unordered map to handle out correspondence in following switch statement
    const static phmap::flat_hash_map<std::string,int> string_to_case{
        {"fasta",1},
        {"fa",1},
        {"fasta.gz",1},
        {"fa.gz",1},
        {"fastq",2},
        {"fq",2},
        {"fastq.gz",2},
        {"fq.gz",2},
        {"bam",2},
        {"gfa",3},
        {"gfa.gz",3},
        {"gfa2",4},
        {"gfa2.gz",4},
        {"vcf",5},
        {"vcf.gz",5}
    };
    
    std::string path = rmFileExt(file), ext = getFileExt("." + file); // variable to handle output path and extension
    
    if(getFileExt(ext) == ".gz") // depending on use input get output format
        gzip = true;
    
    // if the input is not in the unordered map, it means we need to write a new file with the path provided by the user otherwise the output is in the format specified by the user
    if (string_to_case.find(path) == string_to_case.end())
        outFile = true;
    else
        ext = file;
        
    lg.verbose("Writing ouput: " + ext);
    
    if (gzip && outFile) { // if the requested output is gzip compressed and should be outputted to a file
        
        stream = std::make_unique<std::ostream>(zfout.rdbuf()); // then we use the stream for gzip compressed file outputs
        zfout.addHeader();
        
    }else if (!gzip && outFile){ // else if no compression is requested
        
        stream = std::make_unique<std::ostream>(ofs.rdbuf());  // we use the stream regular file outputs
        
    }else{ // else the output is not written to a file
        
        // we close and delete the file
        ofs.close();
        remove(file.c_str());
        
        if (gzip) { // if the output to stdout needs to be compressed we use the appropriate stream
            
            stream = std::make_unique<std::ostream>(zout.rdbuf());
            zout.addHeader();
            
        }else{ // else we use a regular cout stream
            
            std::cout.flush();
            stream = std::make_unique<std::ostream>(std::cout.rdbuf());
        }
    }    
}

OutputStream::~OutputStream() {
    if(this->gzip && this->outFile) // if we wrote to file as gzip, we add the footer and close
        zfout.close();
    else if(gzip && !outFile) // if we streamed as gzip, we add the footer and close
        zout.close();

    if(outFile) // if we wrote to file, we close it
        ofs.close();
}

bool Report::segmentReport(InSequences &inSequences, int &outSequence_flag) { // method to output the summary statistics for each segment
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal poinst
    counter = 1;
    std::vector<InSegment*>* inSegments = inSequences.getInSegments();
    std::vector<unsigned int> circularSegments = inSequences.getCircularSegments();
    
    std::cout<<output("Seq\tHeader\tComment\tTotal segment length\tA\tC\tG\tT\tGC content %\t# soft-masked bases\tIs circular");
    
    if (outSequence_flag)
        std::cout<<output("Sequence\tQuality");
    
    for (InSegment* inSegment : *inSegments) {
        
        std::cout   <<"\n"<<counter++<<"\t"
                    <<inSegment->getSeqHeader()<<"\t"
                    <<inSegment->getSeqComment()<<"\t"
                    <<inSegment->getSegmentLen()<<"\t"
                    <<inSegment->getA()<<"\t"
                    <<inSegment->getC()<<"\t"
                    <<inSegment->getG()<<"\t"
                    <<inSegment->getT()<<"\t"
                    <<inSegment->computeGCcontent()<<"\t"
                    <<inSegment->getLowerCount()<<"\t"
                    <<(inSegment->isCircularSegment(&circularSegments) ? "Y" : "N");

        if (outSequence_flag)
            std::cout<<"\t"<<inSegment->getInSequence()<<"\t"<<inSegment->getInSequenceQuality()<<"\n";
    }
    std::cout<<std::endl;
    return true;
}

bool Report::pathReport(InSequences &inSequences) { // method to output the summary statistics for each path
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal poinst
    counter = 1;
    std::vector<InPath> inPaths = inSequences.getInPaths();
    std::vector<uint32_t> circularPaths = inSequences.getCircularPaths();
    
    std::cout<<output("Seq\tHeader\tComment\tTotal path length\tTotal segment length\tContig N\tA\tC\tG\tT\t# soft-masked bases\tIs circular");
    
    for (InPath& inPath : inPaths) { // loop through all paths
        
        inSequences.walkPath(&inPath);

        std::cout<<"\n"<<counter++<<"\t"
                 <<inPath.getHeader()<<"\t"
                 <<inPath.getComment()<<"\t"
                 <<inPath.getLen()<<"\t"
                 <<inPath.getSegmentLen()<<"\t"
                 <<inPath.getContigN()<<"\t"
                 <<inPath.getA()<<"\t"
                 <<inPath.getC()<<"\t"
                 <<inPath.getG()<<"\t"
                 <<inPath.getT()<<"\t"
                 <<inPath.getLowerCount()<<"\t"
                 <<(inPath.isCircularPath(&circularPaths) ? "Y" : "N");
    }
    std::cout<<std::endl;
    return true;
}

void Report::writeToStream(InSequences &inSequences, std::string file, UserInput &userInput) {
    
    OutputStream outputStream(file);
    
    const static phmap::flat_hash_map<std::string,int> string_to_case{
        {"fasta",1},
        {"fa",1},
        {"fasta.gz",1},
        {"fa.gz",1},
        {"fastq",2},
        {"fq",2},
        {"fastq.gz",2},
        {"fq.gz",2},
        {"gfa",3},
        {"gfa.gz",3},
        {"gfa2",4},
        {"gfa2.gz",4},
        {"vcf",5},
        {"vcf.gz",5}
    };
    
    userInput.stats_flag = outputStream.outFile ? true : false; // since we write to file, let's output the stats
    std::string ext = outputStream.outFile ? getFileExt(outputStream.file) : file; // variable to handle output path and extension
    
    std::unique_ptr<std::ostream> &stream = outputStream.stream;
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) { // this switch allows us to generate the output according to the input request and the unordered map. If the requested output format is not in the map we fall back to the undefined case (0)
            
        case 1: { // fasta[.gz]
            
            std::string pHeader;
            std::string inSeq; // the new sequence being built
            std::vector<InPath> inPaths = inSequences.getInPaths();
            std::vector<InSegment*>* inSegments = inSequences.getInSegments();
            std::vector<InGap>* inGaps = inSequences.getInGaps();
            std::vector<InEdge>* inEdges = inSequences.getEdges();
            std::vector<PathComponent> pathComponents;
            
            unsigned int uId = 0;
            
            for (InPath& inPath : inSequences.getInPaths()) {
                
                if (inPath.getHeader() == "")
                    pHeader = inPath.getpUId();
                else
                    pHeader = inPath.getHeader();
                
                *stream <<">"
                <<pHeader;
                
                if (inPath.getComment() != "") {
                    
                    *stream <<" "
                    <<inPath.getComment();
                }
                
                *stream << std::endl;
                
                pathComponents = inPath.getComponents();
				uint64_t ovlLen = 0;
                
                for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                    
                    uId = component->id;
                    
                    if (component->componentType == SEGMENT) {
                        
                        auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                        if (sId == inSegments->end()) {std::cout<<"Error: cannot find path component. Terminating."<<std::endl; exit(1);} // gives us the segment index
                        
                        if ((*sId)->getInSequencePtr() == NULL) {std::cout<<"Error: Fasta output not possible without segment sequence. Terminating."<<std::endl; exit(0);}
						else if (component->orientation == '+' && (component->end-component->start >= ovlLen || component->end-component->start == 0))
							inSeq += (*sId)->getInSequence(component->start, component->end).substr(ovlLen);
                        else if (component->end-component->start >= ovlLen || component->end-component->start == 0)
                            inSeq += revCom((*sId)->getInSequence(component->start, component->end)).substr(ovlLen);
						else {std::cout<<"Error: overlap longer than component. Terminating."<<std::endl; exit(1);}
                        
                    }else if(component->componentType == EDGE){ // this is just a prototype, need to handle cigar
                        
                        auto edge = find_if(inEdges->begin(), inEdges->end(), [uId](InEdge& obj) {return obj.geteUId() == uId;}); // given a node Uid, find it
                        if (edge == inEdges->end()) {std::cout<<"Error: cannot find path component. Terminating."<<std::endl; exit(0);} // gives us the edge index
                        
						ovlLen = parseCigar(edge->getCigar());
                    }else if(component->componentType == GAP) {
                        
                        auto gap = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                        if (gap == inGaps->end()) {std::cout<<"Error: cannot find path component. Terminating."<<std::endl; exit(1);} // gives us the gap index
                        inSeq += std::string(gap->getDist(component->start, component->end), 'N');
					}else{
						std::cout<<"Error: Unrecognized component type. Terminating."<<std::endl; exit(1);
					}
                }
                if (userInput.splitLength != 0)
                    textWrap(inSeq, *stream, userInput.splitLength); // wrapping text at user-specified line length
                else{
                    *stream<<inSeq<<std::endl;
                    inSeq = "";
                    (*stream).flush();
                }
            }
            break;
        }
        case 2: { // fastq[.gz]
            
            std::string pHeader;
            std::string inSeq, inSeqQual; // the new sequence being built and its quality
            std::vector<InPath> inPaths = inSequences.getInPaths();
            std::vector<InSegment*>* inSegments = inSequences.getInSegments();
            std::vector<InGap>* inGaps = inSequences.getInGaps();
            std::vector<PathComponent> pathComponents;
            
            unsigned int uId = 0, sIdx = 0, gIdx = 0;
            
            for (InPath inPath : inSequences.getInPaths()) {
                
                if (inPath.getHeader() == "")
                    pHeader = inPath.getpUId();
                else
                    pHeader = inPath.getHeader();
                
                *stream <<"@"
                <<pHeader;
                
                if (inPath.getComment() != "")
                    *stream <<" " << inPath.getComment();
                
                *stream <<"\n";
                
                pathComponents = inPath.getComponents();
                
                for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                    
                    uId = component->id;
                    
                    if (component->componentType == SEGMENT) {
                        
                        auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                        
                        if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                        
                        if (component->orientation == '+') {
                            
                            inSeq += (*inSegments)[sIdx]->getInSequence(component->start, component->end);
                            
                        }else{
                            
                            inSeq += revCom((*inSegments)[sIdx]->getInSequence(component->start, component->end));
                            
                        }
                        
                        if ((*inSegments)[sIdx]->getInSequenceQuality() != "") {
                            
                            inSeqQual += (*inSegments)[sIdx]->getInSequenceQuality(component->start, component->end);
                            
                        }else{
                            
                            inSeqQual += std::string((*inSegments)[sIdx]->getInSequence().size(), '!');
                            
                        }
                        
                    }else{
                        
                        auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                        
                        if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the segment index
                        
                        inSeq += std::string((*inGaps)[gIdx].getDist(), 'N');
                        inSeqQual += std::string((*inGaps)[gIdx].getDist(), '!');
                    }
                }
                *stream<<inSeq<<"\n+\n"<<inSeqQual<<"\n";
                
                inSeq = "";
                inSeqQual = "";
            }
            break;
        }
        case 3: { // gfa[.gz] GFA1.2
            
            std::string seqHeader, gHeader, pHeader;
            
            phmap::flat_hash_map<unsigned int, std::string> idsToHeaders = *inSequences.getHash2();
            
            std::vector<InSegment*>* inSegments = inSequences.getInSegments();
            
            *stream<<"H\tVN:Z:1.2\n";
            
            for (InSegment* inSegment : *inSegments) {
                
                seqHeader = inSegment->getSeqHeader();
                
                *stream <<"S\t" // line type
                <<seqHeader<<"\t"; // header
                
                if (userInput.noSequence)
                    *stream <<'*'; // sequence
                else
                    *stream <<inSegment->getInSequence(); // sequence
                
                if (userInput.noSequence)
                    *stream <<"\tLN:i:"<<inSegment->getSegmentLen(); // tags
                
                std::vector<Tag> tags = inSegment->getTags();
                
                for (Tag &tag : tags)
                    *stream <<"\t"<<tag.label<<":"<<tag.type<<":"<<tag.content; // tags
                
                if (inSegment->getInSequenceQuality() != "")
                    *stream <<"\tQL:Z:"<<inSegment->getInSequenceQuality(); // optional quality
                
                *stream<<"\n";
                
            }
            
            for (InEdge inEdge : *(inSequences.getEdges())) {
                
                *stream <<"L\t" // line type
                <<idsToHeaders[inEdge.getsId1()]<<"\t"<<inEdge.getsId1Or()<<"\t"
                <<idsToHeaders[inEdge.getsId2()]<<"\t"<<inEdge.getsId2Or()<<"\t";
                
                *stream <<inEdge.getCigar(); // CIGAR
                
                std::vector<Tag> tags = inEdge.getTags();
                
                for (Tag &tag : tags)
                    *stream <<"\t"<<tag.label<<":"<<tag.type<<":"<<tag.content; // tags
   
                *stream <<"\n";
                
            }
            
            for (InGap inGap : inSequences.getGaps()) {
                
                *stream <<"J\t" // line type
                <<idsToHeaders[inGap.getsId1()]<<"\t"<<inGap.getsId1Or()<<"\t"
                <<idsToHeaders[inGap.getsId2()]<<"\t"<<inGap.getsId2Or()<<"\t";
                
                if (inGap.getDist() != 0)
                    *stream <<inGap.getDist(); // gap size
                else
                    *stream <<"*"; // gap size

                std::vector<Tag> tags = inGap.getTags();
                
                for (Tag &tag : tags)
                    *stream <<"\t"<<tag.label<<":"<<tag.type<<":"<<tag.content; // tags
                
                *stream <<"\n";
            }
            std::vector<PathComponent> pathComponents;
            
            for (InPath inPath : inSequences.getInPaths()) {
                
                if (inPath.getHeader() == "")
                    pHeader = inPath.getpUId();
                else
                    pHeader = inPath.getHeader();
    
                *stream <<"P\t" // line type
                <<pHeader<<"\t"; // id
                
                pathComponents = inPath.getComponents();
                
                for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                    
                    if(component->orientation != '0') {
                        
                        *stream << idsToHeaders[component->id];
                        
                        if(component->start != 0 || component->end != 0)
                            *stream << "(" << std::to_string(component->start) << ":" << std::to_string(component->end) << ")";
                        
                        *stream << component->orientation;
                        
                    }else
                        *stream <<(component->componentType == EDGE ? ',' : ';');
                }
                
                *stream << "\t*"; // cigar
                
                if (inPath.getComment() != "")
                    *stream <<"\tCM:Z:"<<inPath.getComment();
                *stream <<"\n";
            }
            break;
        }
        case 4: { // gfa2[.gz] GFA2
            
            std::string seqHeader, gHeader, pHeader;
            phmap::flat_hash_map<unsigned int, std::string> idsToHeaders = *inSequences.getHash2();
            std::vector<InSegment*>* inSegments = inSequences.getInSegments();
            
            *stream<<"H\tVN:Z:2.0\n";
            
            for (InSegment* inSegment : *inSegments) {
                
                seqHeader = inSegment->getSeqHeader();
                
                *stream <<"S\t" // line type
                <<seqHeader<<"\t" // header
                <<inSegment->getSegmentLen()<<"\t" // seq length
                <<inSegment->getInSequence(); // sequence
                
                std::vector<Tag> tags = inSegment->getTags();
                
                for (Tag &tag : tags)
                    *stream <<"\t"<<tag.label<<":"<<tag.type<<":"<<tag.content; // tags
                
                if (inSegment->getInSequenceQuality() != "")
                    *stream <<"\tQL:Z:"<<inSegment->getInSequenceQuality(); // optional quality
                
                *stream<<"\n";
            }
            
            for (InEdge inEdge : *(inSequences.getEdges())) {
                
                *stream <<"E\t" // line type
                <<inEdge.geteHeader()<<"\t"
                <<idsToHeaders[inEdge.getsId1()]<<"\t"<<inEdge.getsId1Or()<<"\t" // sUid1:sid1:ref
                <<idsToHeaders[inEdge.getsId2()]<<"\t"<<inEdge.getsId2Or()<<"\t"; // sUid2:sid2:ref
                
                *stream <<inEdge.getCigar(); // CIGAR
                
                std::vector<Tag> tags = inEdge.getTags();
                
                for (Tag &tag : tags)
                    *stream <<"\t"<<tag.label<<":"<<tag.type<<":"<<tag.content; // tags
                
                *stream <<"\n";
            }
            
            for (InGap inGap : inSequences.getGaps()) {
                
                if (inGap.getgHeader() == "")
                    gHeader = inGap.getuId();
                else
                    gHeader = inGap.getgHeader();
                
                *stream <<"G\t" // line type
                <<gHeader<<"\t" // id
                <<idsToHeaders[inGap.getsId1()]<<inGap.getsId1Or()<<"\t" // sUid1:sid1:ref
                <<idsToHeaders[inGap.getsId2()]<<inGap.getsId2Or()<<"\t" // sUid2:sid2:ref
                <<inGap.getDist(); // size
                
                std::vector<Tag> tags = inGap.getTags();
                
                for (Tag &tag : tags)
                    *stream <<"\t"<<tag.label<<":"<<tag.type<<":"<<tag.content; // tags
                
                *stream <<"\n";
            }
            std::vector<PathComponent> pathComponents;
            
            for (InPath inPath : inSequences.getInPaths()) {
                
                if (inPath.getHeader() == "")
                    pHeader = inPath.getpUId();
                else
                    pHeader = inPath.getHeader();
                    
                *stream <<"O\t" // line type
                <<pHeader<<"\t"; // id
                
                pathComponents = inPath.getComponents();
                
                for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                    
                    *stream << idsToHeaders[component->id];
                    
                    if(component->start != 0 || component->end != 0)
                        *stream << "(" << std::to_string(component->start) << ":" << std::to_string(component->end) << ")";
                    
                    if(component->orientation != '0')
                        *stream << component->orientation;
                    
                    if (component != std::prev(pathComponents.end()))
                        *stream <<" "; // space
                }
                
                if (inPath.getComment() != "") {
                    
                    *stream <<"\tCM:Z:"
                    <<inPath.getComment();
                }
                *stream <<"\n";
            }
            break;
        }
        case 5: { // .vcf[.gz]
            
            std::string pHeader;
            std::vector<InSegment*> *inSegments = inSequences.getInSegments();
            std::vector<InGap> *inGaps = inSequences.getInGaps();
            
            *stream<<"##fileformat=VCFv4.2\n";
            *stream<<"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
            *stream<<"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
            *stream<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
            
            //            for (InPath& inPath : inSequences.getInPaths())
            //                *stream<<"contigs\n";
            
            for (InPath& inPath : inSequences.getInPaths()) {
                
                if (inPath.getHeader() == "")
                    pHeader = inPath.getpUId();
                else
                    pHeader = inPath.getHeader();
                
                unsigned int cUId = 0, gapLen = 0;
                std::vector<PathComponent> pathComponents = inPath.getComponents();
                uint64_t absPos = 1; // vcf is 1-based
                
                for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                    
                    cUId = component->id;
                    
                    if (component->componentType == SEGMENT) {
                        
                        auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                        std::string* ref = (*inSegment)->getInSequencePtr();
                        
                        if (component->orientation == '+') {
                            
                            std::vector<std::deque<DBGpath>> variants = (*inSegment)->getVariants();
                            
                            for (std::deque<DBGpath> DBGpaths : variants) {
                                
                                uint64_t pos = DBGpaths[0].pos;
                                bool indel = false;
                                
                                if (std::any_of(DBGpaths.begin(), DBGpaths.end(), [](DBGpath& obj){return obj.type == DEL || obj.type == INS;})) // indel
                                    indel = true;
                                
                                uint8_t offset = 0;
                                if (indel)
                                    offset = 1;
                                
                                std::string refSeq;
                                for (uint16_t b = 0; b < DBGpaths[0].refLen; ++b)
                                    refSeq.push_back((*ref)[pos-offset+b]);
                                
                                if (indel)
                                    *stream<<pHeader<<"\t"<<absPos+pos-1<<"\t.\t"<<refSeq<<(*ref)[pos+DBGpaths[0].refLen-1]<<"\t";
                                else
                                    *stream<<pHeader<<"\t"<<absPos+pos<<"\t.\t"<<refSeq<<"\t";
                                
                                double qual = 0;
                                bool first = true;
                                for (const DBGpath& variant : DBGpaths) {
                                    
                                    if (first) {first = false;} else {*stream<<",";}
                                    
                                    if(variant.type == SNV || variant.type == COM) {
                                        if (indel)
                                            *stream<<(*ref)[pos-1]<<variant.sequence;
                                        else
                                            *stream<<variant.sequence;
                                        
                                    }else if(variant.type == DEL) {
                                        *stream<<(*ref)[pos-1]<<variant.sequence<<(*ref)[pos];
                                    }else if(variant.type == INS) {
                                        *stream<<(*ref)[pos-1];
                                    }
                                    
                                    qual += variant.score;
                                    
                                }
                                *stream<<"\t"<<round(qual/DBGpaths.size())<<"\tPASS\t.\tGT:GQ\t1/1:"
                                <<DBGpaths[0].score
                                <<"\n";
                            }
                            absPos += (*inSegment)->getSegmentLen();
                            
                        }else{} // GFA not handled yet
                    }else if (component->componentType == GAP){
                        
                        auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                        gapLen = inGap->getDist(component->start - component->end);
                        absPos += gapLen;
                        
                    }else{} // need to handle edges, cigars etc
                }
            }
            break;
        }
        case 0: { // undefined case
            std::cout<<"Unrecognized output format: "<<outputStream.file;
            break;
        }
    }
}

bool Report::outCoord(InSequences &inSequences, char bedOutType, bool sizeOnly) { // method to output the coordinates of each feature
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal poinst

    counter = 0;
    
    std::string seqHeader;
    std::vector<unsigned int> seqBoundaries;
    
    unsigned int uId = 0, sIdx = 0, gIdx = 0;
    uint64_t pos = 0;

    std::string pHeader;
    std::vector<InPath> inPaths = inSequences.getInPaths();
    std::vector<InSegment*>* inSegments = inSequences.getInSegments();
    std::vector<InGap>* inGaps = inSequences.getInGaps();
    std::vector<PathComponent> pathComponents;

    for (InPath inPath : inSequences.getInPaths()) {
        
        if (inPath.getHeader() == "") {
            
            pHeader = inPath.getpUId();
            
        }else{
            
            pHeader = inPath.getHeader();
            
        }

        pathComponents = inPath.getComponents();

        switch (bedOutType) {
                
            case 's': { // scaffolds
                
                std::cout<<pHeader<<"\t";
                
                if (!sizeOnly) {
                
                    std::cout<<pos<<"\t";
                    
                }
                
                for (auto &component : pathComponents) {
                    
                    uId = component.id;
                    
                    if (component.componentType == SEGMENT) {
                    
                        auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                        
                        if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                        
                        pos += (*inSegments)[sIdx]->getInSequence().size();
                        
                    }else{
                        
                        auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                        
                        if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the gap index
                        
                        pos += (*inGaps)[gIdx].getDist();
                        
                    }
                    
                    
                }
                
                std::cout<<pos<<"\n";
                pos = 0;
                
                break;
            }

            case 'c': { // contigs
                
                for (auto &component : pathComponents) {
                    
                    uId = component.id;
                    
                    if (component.componentType == SEGMENT) {
                    
                        auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                        
                        if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                        
                        std::cout<<(sizeOnly ? (*inSegments)[sIdx]->getSeqHeader() : pHeader)<<"\t";
                        
                        if (!sizeOnly) {
                        
                            std::cout<<pos<<"\t";
                            
                        }
                        
                        if (sizeOnly) {
                            
                            pos = (*inSegments)[sIdx]->getInSequence().size();
                            
                        }else{
                            
                            pos += (*inSegments)[sIdx]->getInSequence().size();
                        
                        }
                        
                        std::cout<<pos<<"\n";
                        
                        if (sizeOnly)
                            pos = 0;
                        
                    }else{
                        
                        auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                        
                        if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the gap index
                        
                        pos += (*inGaps)[gIdx].getDist();
                        
                    }
                }
                
                pos = 0;
                
                break;
            }
            
            case 'g': { // gaps
                
                for (auto &component : pathComponents) {
                    
                    uId = component.id;
                    
                    if (component.componentType == SEGMENT) {
                    
                        auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                        
                        if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                        
                        pos += (*inSegments)[sIdx]->getInSequence().size();
                        
                    }else{
                        
                        auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                        
                        if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the gap index
                        
                        std::cout<<(sizeOnly ? (*inGaps)[gIdx].getgHeader() : pHeader)<<"\t";
                        
                        if (!sizeOnly) {
                        
                            std::cout<<pos<<"\t";
                            
                        }
                        
                        if (sizeOnly) {
                            
                            pos = (*inGaps)[gIdx].getDist();
                            
                        }else{
                            
                            pos += (*inGaps)[gIdx].getDist();
                        
                        }
                            
                        std::cout<<pos<<"\n";
                        
                    }
                    
                }
                
                pos = 0;
                
                break;
            }

            case 'h': { // homopolymer runs

                for (auto &component : pathComponents) {
                    
                    uId = component.id;
                    
                    if (component.componentType == SEGMENT) {

                        auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                        
                        if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index

                        InSegment* seg = inSegments->at(sIdx);

                        std::vector<std::pair<unsigned int, unsigned int>> bedCoords;
                        homopolymerBedCoords(seg->inSequence, bedCoords, 1);

                        for(const auto &pair : bedCoords) {
                            std::cout << pHeader << "\t" << pair.first+pos << "\t" << pair.second+pos << std::endl;
                        }

                        pos += seg->inSequence->size();

                    }else{
                        
                        auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                        
                        if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the gap index
                        
                        pos += (*inGaps)[gIdx].getDist();
                        
                    }
                    
                }
                
                pos = 0;

                break;
            }
            
            // default includes 'a'
            default: { // both contigs and gaps in .agp format
                
                for (auto &component : pathComponents) {
                    
                    uId = component.id;
                    
                    if (component.componentType == SEGMENT) {
                    
                        auto sId = find_if(inSegments->begin(), inSegments->end(), [uId](InSegment* obj) {return obj->getuId() == uId;}); // given a node Uid, find it
                        
                        if (sId != inSegments->end()) {sIdx = std::distance(inSegments->begin(), sId);} // gives us the segment index
                        
                        std::cout<<pHeader<<"\t"<<pos+1;
                        
                        pos += (*inSegments)[sIdx]->getInSequence().size();
                        
                        counter++;
                        
                        std::cout<<"\t"<<pos<<"\t"<<counter<<"\tW\t"<<(*inSegments)[sIdx]->getSeqHeader()<<"\t1\t"<<(*inSegments)[sIdx]->getInSequence().size()<<"\t"<<component.orientation<<"\n";
                        
                    }else{
                        
                        auto gId = find_if(inGaps->begin(), inGaps->end(), [uId](InGap& obj) {return obj.getuId() == uId;}); // given a node Uid, find it
                        
                        if (gId != inGaps->end()) {gIdx = std::distance(inGaps->begin(), gId);} // gives us the gap index
                        
                        std::cout<<pHeader<<"\t"<<pos+1;
                        
                        pos += (*inGaps)[gIdx].getDist();
                        
                        counter++;
                        
                        std::cout<<"\t"<<pos<<"\t"<<counter<<"\tN\t"<<(*inGaps)[gIdx].getDist()<<"\t"<<(*inGaps)[gIdx].getgHeader()<<"\tyes\n";
                        
                    }
                    
                }
                
                pos = 0;
                counter = 0;
            }
        }
    }

    return true;
}

bool Report::reportStats(InSequences &inSequences, uint64_t gSize, int outBubbles_flag) { // method to output all summary statistics for the entire sequence set
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal poinst

    if (!tabular_flag) {
    
        std::cout<<output("+++Assembly summary+++")<<"\n";
    
    }
    
    if (gSize > 0) {
    
        std::cout<<output("Expected genome size")<<gSize<<"\n";
    
    }
    
    std::cout<<output("# scaffolds")<<inSequences.getScaffN()<<"\n";
    std::cout<<output("Total scaffold length")<<inSequences.getTotScaffLen()<<"\n";
    std::cout<<output("Average scaffold length") << gfa_round(inSequences.computeAvgScaffLen()) << "\n";
    inSequences.evalNstars('s', gSize); // scaffold N* statistics
    std::cout<<output("Scaffold N50")<<inSequences.getScaffN50()<<"\n";
    inSequences.evalAuN('s', gSize); // scaffold auN
    std::cout<<output("Scaffold auN") << gfa_round(inSequences.getScaffauN()) << "\n";
    std::cout<<output("Scaffold L50")<<inSequences.getScaffL50()<<"\n";
    
    if (gSize > 0) {
        
        std::cout<<output("Scaffold NG50")<<inSequences.getScaffNG50()<<"\n";
        std::cout<<output("Scaffold auNG") << gfa_round(inSequences.getScaffauNG()) << "\n";
        std::cout<<output("Scaffold LG50")<<inSequences.getScaffLG50()<<"\n";
        
    }
    std::cout<<output("Largest scaffold")<<inSequences.getLargestScaffold()<<"\n";
    std::cout<<output("Smallest scaffold")<<inSequences.getSmallestScaffold()<<"\n";
    
    std::cout<<output("# contigs")<<inSequences.getTotContigN()<<"\n";
    std::cout<<output("Total contig length")<<inSequences.getTotContigLen()<<"\n";
    std::cout<<output("Average contig length") << gfa_round(inSequences.computeAvgContigLen()) << "\n";
    inSequences.evalNstars('c', gSize); // contig N* statistics
    std::cout<<output("Contig N50")<<inSequences.getContigN50()<<"\n";
    inSequences.evalAuN('c', gSize); // contig auN
    std::cout<<output("Contig auN") << gfa_round(inSequences.getContigauN()) << "\n";
    std::cout<<output("Contig L50")<<inSequences.getContigL50()<<"\n";
    
    if (gSize > 0) {
        
        std::cout<<output("Contig NG50")<<inSequences.getContigNG50()<<"\n";
        std::cout<<output("Contig auNG") << gfa_round(inSequences.getContigauNG()) << "\n";
        std::cout<<output("Contig LG50")<<inSequences.getContigLG50()<<"\n";
        
    }
    std::cout<<output("Largest contig")<<inSequences.getLargestContig()<<"\n";
    std::cout<<output("Smallest contig")<<inSequences.getSmallestContig()<<"\n";
    
    std::cout<<output("# gaps in scaffolds")<<inSequences.getGapNScaffold()<<"\n";
    std::cout<<output("Total gap length in scaffolds")<<inSequences.getTotGapLen()<<"\n";
    std::cout<<output("Average gap length in scaffolds") << gfa_round(inSequences.computeAverageGapLen()) << "\n";
    inSequences.evalNstars('g'); // gap N* statistics
    std::cout<<output("Gap N50 in scaffolds")<<inSequences.getGapN50()<<"\n";
    inSequences.evalAuN('g'); // gap auN
    std::cout<<output("Gap auN in scaffolds") << gfa_round(inSequences.getGapauN()) << "\n";
    std::cout<<output("Gap L50 in scaffolds")<<inSequences.getGapL50()<<"\n";
    std::cout<<output("Largest gap in scaffolds")<<inSequences.getLargestGap()<<"\n";
    std::cout<<output("Smallest gap in scaffolds")<<inSequences.getSmallestGap()<<"\n";
    
    std::cout<<output("Base composition (A:C:G:T)");
    std::cout << inSequences.getTotA() << ":"
              << inSequences.getTotC() << ":"
              << inSequences.getTotG() << ":"
              << inSequences.getTotT() << "\n";
    std::cout<<output("GC content %") << gfa_round(inSequences.computeGCcontent()) << "\n";
    std::cout<<output("# soft-masked bases")<<inSequences.getTotLowerCount()<<"\n";
    
    // graph statistics
    std::cout<<output("# segments")<<inSequences.getSegmentN()<<"\n";
    std::cout<<output("Total segment length")<<inSequences.getTotSegmentLen()<<"\n";
    std::cout<<output("Average segment length") << gfa_round(inSequences.computeAvgSegmentLen()) << "\n";
    
    std::cout<<output("# gaps")<<inSequences.getGapN()<<"\n";
    std::cout<<output("# paths")<<inSequences.getPathN()<<"\n";
    
    counter = 0;
    unsigned int connectedComponents = 0;
    
    unsigned int edgeN = inSequences.getEdgeN();

    if (edgeN > 0) {
        
        std::cout<<output("# edges")<<edgeN<<"\n";
        std::cout<<output("Average degree")<<(double)inSequences.getEdgeN()/inSequences.getSegmentN()<<"\n";
    
        inSequences.buildEdgeGraph();

        lg.verbose("Graph DFS");
        
        std::vector<InSegment*>* inSegments = inSequences.getInSegments();
        std::vector<unsigned int> componentLengths;
        unsigned int componentLength = 0;
        
        for (InSegment* inSegment : *inSegments) { // loop through all nodes
            
            if (!inSequences.getVisited(inSegment->getuId()) && !inSequences.getDeleted(inSegment->getuId())) { // check if the node was already visited
                
                inSequences.dfsEdges(inSegment->getuId(), &componentLength); // if not, visit all connected components recursively
                ++connectedComponents;
                componentLengths.push_back(componentLength);
                componentLength = 0;
            }
        }
        sort(componentLengths.begin(), componentLengths.end(), std::greater<unsigned int>());
        std::cout<<output("# connected components")<<connectedComponents-inSequences.getDisconnectedComponents()<<"\n";
        std::cout<<output("Largest connected component length")<<componentLengths[0]<<"\n";
        std::cout<<output("# dead ends")<<inSequences.getDeadEnds()<<"\n";
        std::cout<<output("# disconnected components")<<inSequences.getDisconnectedComponents()<<"\n";
        std::cout<<output("Total length disconnected components")<<inSequences.getLengthDisconnectedComponents()<<"\n";
        std::cout<<output("# separated components")<<connectedComponents<<"\n";
        
        inSequences.findBubbles();
        
        std::cout<<output("# bubbles")<<inSequences.getBubbles()->size()<<"\n";
        
        if (outBubbles_flag) {
            
            phmap::flat_hash_map<unsigned int, std::string> idsToHeaders = *inSequences.getHash2();
            
            for (Bubble bubble : *inSequences.getBubbles()) { // loop through all nodes
                
                std::cout<<idsToHeaders[bubble.id0]<<"\t"
                         <<idsToHeaders[bubble.id1]<<"\t"
                         <<idsToHeaders[bubble.id2]<<"\t"
                         <<idsToHeaders[bubble.id3]<<"\n";
                
            }
        
        }
        
        std::cout<<output("# circular segments")<<inSequences.getCircularSegments().size()<<"\n";
        std::cout<<output("# circular paths")<<inSequences.getCircularPaths().size()<<"\n";
            
    }

    return true;
    
}

bool Report::nstarReport(InSequences &inSequences, uint64_t gSize) { // method to generate all N** reports
    std::cout << std::fixed; // disables scientific notation
    std::cout << std::setprecision(2); // 2 decimal poinst

    int pos = 1;
    std::vector <uint64_t> scaffNstars = inSequences.getScaffNstars();
    for (unsigned int val : scaffNstars) {
        std::cout<<output("Scaffold N"+std::to_string(pos*10))<<val<<"\n";
        pos++;
    }
    
    pos = 1;
    std::vector <unsigned int> scaffLstars = inSequences.getScaffLstars();
    for (unsigned int val : scaffLstars) {
        std::cout<<output("Scaffold L"+std::to_string(pos*10))<<val<<"\n";
        pos++;
    }
    
    if (gSize > 0) {
        
        pos = 1;
        std::vector <uint64_t> scaffNGstars = inSequences.getScaffNGstars();
        for (unsigned int val : scaffNGstars) {
            std::cout<<output("Scaffold NG"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> scaffLGstars = inSequences.getScaffLGstars();
        for (unsigned int val : scaffLGstars) {
            std::cout<<output("Scaffold LG"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
    }
    
    pos = 1;
    std::vector <uint64_t> contigNstars = inSequences.getContigNstars();
    for (unsigned int val : contigNstars) {
        std::cout<<output("Contig N"+std::to_string(pos*10))<<val<<"\n";
        pos++;
    }
    
    pos = 1;
    std::vector <unsigned int> contigLstars = inSequences.getContigLstars();
    for (unsigned int val : contigLstars) {
        std::cout<<output("Contig L"+std::to_string(pos*10))<<val<<"\n";
        pos++;
    }
    
    if (gSize > 0) {
        
        pos = 1;
        std::vector <uint64_t> contigNGstars = inSequences.getContigNGstars();
        for (unsigned int val : contigNGstars) {
            std::cout<<output("Contig NG"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> contigLGstars = inSequences.getContigLGstars();
        for (unsigned int val : contigLGstars) {
            std::cout<<output("Contig LG"+std::to_string(pos*10))<<val<<"\n";
            pos++;
        }
        
    }
    
    pos = 1;
    std::vector <uint64_t> gapNstars = inSequences.getGapNstars();
    for (unsigned int val : gapNstars) {
        std::cout<<output("Gap N"+std::to_string(pos*10))<<val<<"\n";
        pos++;
    }
    
    pos = 1;
    std::vector <unsigned int> gapLstars = inSequences.getGapLstars();
    for (unsigned int val : gapLstars) {
        std::cout<<output("Gap L"+std::to_string(pos*10))<<val<<"\n";
        pos++;
    }
    
    return true;
    
}
