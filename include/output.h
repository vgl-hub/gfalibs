//
//  gfastats-output.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef OUTPUT_H
#define OUTPUT_H

//classes
class Report {

private:
    unsigned int counter = 0;
    
public:
    bool seqReport(InSequences &inSequences, int &outSequence_flag);
    
    bool outFile(InSequences &inSequences, std::string file, UserInput &userInput, int splitLength = 0);
    
    bool outCoord(InSequences &inSequences, char bedOutType, bool sizeOnly = false);
    
    bool reportStats(InSequences &inSequences, uint64_t gSize, int outBubbles_flag);
    
    bool nstarReport(InSequences &inSequences, uint64_t gSize);
    
};

#endif /* OUTPUT_H */
