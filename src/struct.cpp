#include <cstdint>
#include <string>
#include <vector>

#include "bed.h"
#include "struct.h"

std::string UserInput::file(char type, uint16_t num) {
    
    std::string filename;
    
    switch (type) {
        case 'f':
            filename = inSequence;
            break;
        case 'r':
            filename = inFiles[num];
            break;
        case 'i':
            filename = inBedInclude;
            break;
        case 'e':
            filename = inBedExclude;
            break;
        case 'a':
            filename = inAgp;
            break;
        case 'k':
            filename = inSak;
            break;
        case 'g':
            filename = inAlign;
            break;
    }
    
    return filename;
    
}

Sequences::~Sequences() {
    for (Sequence* sequence : sequences)
        delete sequence;
}

Sequence::~Sequence() {
    if (sequence != NULL)
        delete sequence;
    if (sequenceQuality != NULL)
        delete sequenceQuality;
}

void Sequence::deleteSequence() {
    if (sequence != NULL) {
        delete sequence;
        sequence = NULL;
    }
    if (sequenceQuality != NULL) {
        delete sequenceQuality;
        sequenceQuality = NULL;
    }
}

bool Edge::operator==(const Edge& e) const {
    return orientation0 == e.orientation0 && orientation1 == e.orientation1 && id == e.id;
}
