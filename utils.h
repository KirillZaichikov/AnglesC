#include "consts.h"
#ifndef CHAOS_UTILS_CUH
#define CHAOS_UTILS_CUH

struct LyapunovExponentsAccumulator{
    double lyapunovExp[MAX_DIMENSION];
    unsigned samplesNum;
};

void accumulateExponents(LyapunovExponentsAccumulator* la, double* lyapunovExp, unsigned expNum) {
    for (int i=0; i<expNum; i++) {
        la->lyapunovExp[i] += lyapunovExp[i];
    }
    ++la->samplesNum;
}

void commitExponentAccumulator(LyapunovExponentsAccumulator* laIn, double* lyapunovExpOut, unsigned expNum) {
    for (int i=0; i<expNum; i++) {
        laIn->lyapunovExp[i] = laIn->lyapunovExp[i] / laIn->samplesNum;
    }
}

#endif