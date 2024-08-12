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

struct MemoryForIntegrator{
    double k1[MAX_DIMENSION];
    double k2[MAX_DIMENSION];
    double k3[MAX_DIMENSION];
    double k4[MAX_DIMENSION];
    double k5[MAX_DIMENSION];
    double k6[MAX_DIMENSION];
    double k7[MAX_DIMENSION];
    double k8[MAX_DIMENSION];
    double projSum[MAX_DIMENSION];
    double arg[MAX_DIMENSION];
};

struct MemoryForMatrix{
    double eigenValues[MAX_DIMENSION];
    double A[MAX_DIMENSION*MAX_DIMENSION];
    double R[MAX_DIMENSION*MAX_DIMENSION];
    double Q[MAX_DIMENSION*MAX_DIMENSION];
    double dotProductsMatrix[MAX_DIMENSION * MAX_DIMENSION];
    double dotProductsMatrixT[MAX_DIMENSION * MAX_DIMENSION];
};

void commitExponentAccumulator(LyapunovExponentsAccumulator* laIn, double* lyapunovExpOut, unsigned expNum) {
    for (int i=0; i<expNum; i++) {
        laIn->lyapunovExp[i] = laIn->lyapunovExp[i] / laIn->samplesNum;
    }
}

#endif