
#ifndef _BITWISEDISTANCE_H_
#define _BITWISEDISTANCE_H_

#include <intrin.h>

class BitwiseDistance {

public:

    int countOne(unsigned long bitString) {
        return __popcnt64(bitString);
    }
    int getHammingDistance(unsigned long a, unsigned long b) {
        return countOne(a^b);
    }

};


#endif
