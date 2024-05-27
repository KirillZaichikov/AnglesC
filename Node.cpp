#include "consts.h"
#include "integrator.cpp"
#include "linal.cpp"
#include <cstring>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>


//work
void cudaAngleNodeWarm(double *main, double *sub, double *essCur, const int dimension,
                       void(*diffFunc)(const double*, double*, const double*), void(*diffFuncVar)(const double*, double*, const double*, const double*),
                       void(*diffFuncVarTransRev)(const double*, double*, const double*, const double*),
                       const double *params, double eps, int expsNum, double step, double timeSkip, double timeSkipClP, const int eSSdim) {
    double projSum[MAX_DIMENSION], arg[MAX_DIMENSION];
    double k1[MAX_DIMENSION], k2[MAX_DIMENSION], k3[MAX_DIMENSION], k4[MAX_DIMENSION], k5[MAX_DIMENSION], k6[MAX_DIMENSION], k7[MAX_DIMENSION], k8[MAX_DIMENSION];
    memset(k1, 0., MAX_DIMENSION * sizeof(double));
    memset(k2, 0., MAX_DIMENSION * sizeof(double));
    memset(k3, 0., MAX_DIMENSION * sizeof(double));
    memset(k4, 0., MAX_DIMENSION * sizeof(double));
    memset(k5, 0., MAX_DIMENSION * sizeof(double));
    memset(k6, 0., MAX_DIMENSION * sizeof(double));
    memset(k7, 0., MAX_DIMENSION * sizeof(double));
    memset(k8, 0., MAX_DIMENSION * sizeof(double));
    memset(arg, 0., MAX_DIMENSION * sizeof(double));
    memset(projSum, 0., MAX_DIMENSION * sizeof(double));
    double Twarm = timeSkip / step;
    double Tnorm = timeSkipClP / step;
    printf("TIME FOR SKIP: %lf\n", Twarm);
    printf("TIME FOR NORW: %lf\n", Tnorm);
    for (int t = 0; t < Twarm; t += 1) // Skip points to reach the attractor{
        dverkStep(main, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8); // work

    memset(sub, 0., expsNum * dimension * sizeof(double));
    memset(essCur, 0., expsNum * dimension * sizeof(double));


    for (int i = 0; i < expsNum; i++) {
        sub[i * dimension + i] += eps;
    }
    for (int i = 0; i < expsNum; i++) {
        essCur[i * dimension + i] += eps;
    } //work

    for (int t = 0; t < Tnorm; t += 1) {
        for (int i = 0; i < expsNum; i++)
            dverkStepVarMat(&sub[i * dimension], dimension, diffFuncVar,
                            params, step, arg, main, k1, k2, k3, k4, k5, k6, k7, k8);
        for (int i = 0; i < expsNum; i++)
            dverkStepVarMat(&essCur[i * dimension], dimension, diffFuncVarTransRev, params, step, arg, main, k1, k2, k3, k4, k5, k6, k7, k8);
        dverkStep(main, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        ortVecs(sub, dimension, expsNum, projSum);
        ortVecs(essCur, dimension, expsNum, projSum);
        normalizeVecs(sub, dimension, expsNum, eps);
        normalizeVecs(essCur, dimension, expsNum, eps);
    } //work
}

int cudaAngleNodeCalculate(double *mainTrajectory, double *slaveTrajectories, double *essCur, const int dimension,
                        void(*diffFunc)(const double*, double*,const double*), void(*diffFuncVar)(const double*, double*,const double*,const double*),
                        void(*diffFuncVarRev)(const double*, double*,const double*,const double*),
                        void(*diffFuncVarTransRev)(const double*, double*,const double*,const double*),
                        const double *params, double eps, double *lyapExp, double *lyapExpBw, double *lyapExpT, double *phi,
                        int expsNum, double step, double timeSkip, double timeSkipClP, double calcTime,
                        double addSkipSteps, double *traj_dots, double *ncu_s, const int eSSdim) {
    double projSum[MAX_DIMENSION], arg[MAX_DIMENSION];

    double main_afterForward[MAX_DIMENSION]; // to calc main traj after fw time is completed
    double slaveTrajBw[MAX_DIMENSION * MAX_DIMENSION]; // exp_num * max_dimension
    double dotProductsMatrix[MAX_DIMENSION * MAX_DIMENSION];
    double dotProductsMatrixT[MAX_DIMENSION * MAX_DIMENSION];
    double A[MAX_DIMENSION * MAX_DIMENSION];
    double R[MAX_DIMENSION * MAX_DIMENSION];
    double Q[MAX_DIMENSION * MAX_DIMENSION];
    double eigenValues[MAX_DIMENSION];
    memset(dotProductsMatrix, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(dotProductsMatrixT, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(A, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(R, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(Q, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(eigenValues, 0., MAX_DIMENSION * sizeof(double));

    memset(main_afterForward, 0., MAX_DIMENSION * sizeof(double));
    memset(lyapExp, 0., MAX_DIMENSION * sizeof(double)); //В комплексе поломаный memset
    memset(lyapExpBw, 0., MAX_DIMENSION * sizeof(double));
    memset(lyapExpT, 0., MAX_DIMENSION * sizeof(double));

    int stepCount = 0;
    double k1[MAX_DIMENSION], k2[MAX_DIMENSION], k3[MAX_DIMENSION], k4[MAX_DIMENSION], k5[MAX_DIMENSION], k6[MAX_DIMENSION], k7[MAX_DIMENSION], k8[MAX_DIMENSION];
    memset(k1, 0., MAX_DIMENSION * sizeof(double));
    memset(k2, 0., MAX_DIMENSION * sizeof(double));
    memset(k3, 0., MAX_DIMENSION * sizeof(double));
    memset(k4, 0., MAX_DIMENSION * sizeof(double));
    memset(k5, 0., MAX_DIMENSION * sizeof(double));
    memset(k6, 0., MAX_DIMENSION * sizeof(double));
    memset(k7, 0., MAX_DIMENSION * sizeof(double));
    memset(k8, 0., MAX_DIMENSION * sizeof(double));
    memset(arg, 0., MAX_DIMENSION * sizeof(double));
    memset(projSum, 0., MAX_DIMENSION * sizeof(double));
    double T = calcTime / step; //надо поменять в комплексе на дабл
    double Tnorm = timeSkipClP / step;
    for (int t = 0; t < T; t += 1) {
        for (int i = 0; i < expsNum; i++)
            dverkStepVarMat(&slaveTrajectories[i * dimension], dimension, diffFuncVar, params, step, arg,
                mainTrajectory, k1, k2, k3, k4, k5, k6, k7, k8);
        for (int i = 0; i < expsNum; i++)
            dverkStepVarMat(&essCur[i * dimension], dimension, diffFuncVarTransRev, params, step, arg,
                mainTrajectory, k1, k2, k3, k4, k5, k6, k7, k8);

        ortVecs(slaveTrajectories, dimension, expsNum, projSum);
        ortVecs(essCur, dimension, expsNum, projSum);

        // for (int i = 0; i < dimension; i++)
        //     //printf("%lf ", (traj_dots + stepCount * dimension)[i]);
        //     printf("%lf ", (mainTrajectory)[i]);
        // std::cout << std::endl;
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);

        for (int i = 0; i < expsNum; i++)
            lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / eps);
        for (int i = 0; i < expsNum; i++)
            lyapExpT[i] += log(vecNorm(&essCur[i * dimension], dimension) / eps);

        normalizeVecs(slaveTrajectories, dimension, expsNum, eps);
        normalizeVecs(essCur, dimension, expsNum, eps);
        /*
        if (t % 100 == 0){
            for (int i = 0; i < eSSdim; i++){printf("Ess:");
                for (int j = 0; j < dimension; j++)
                    printf("%lf ", essCur[i*dimension + j]);
                printf("\n");}
            for (int i = 0; i < expsNum; i++){printf("Srr:");
                for (int j = 0; j < dimension; j++)
                    printf("%lf ", slaveTrajectories[i*dimension + j]);
                printf("\n");}
        }*/ //ПОДОЗРИТЕЛЬНО ПОХОЖИ ВЕКТОРА ИЗ SLAVE И ESSCUR
        memcpy(ncu_s + stepCount * dimension * eSSdim, essCur, sizeof(double) * dimension * eSSdim);
        memcpy(traj_dots + stepCount * dimension, mainTrajectory, sizeof(double) * dimension);
        stepCount++;
    } // вроде все ок

    for (int i = 0; i < expsNum; i++)
        lyapExp[i] = lyapExp[i] / calcTime;
    for (int i = 0; i < expsNum; i++)
        lyapExpT[i] = lyapExpT[i] / calcTime;

    memcpy(main_afterForward, mainTrajectory, sizeof(double) * dimension);
    for (int t = 0; t < Tnorm; t += 1) {
        dverkStep(main_afterForward, dimension, diffFunc,
            params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        memcpy(traj_dots + stepCount * dimension, main_afterForward, sizeof(double) * dimension);
        stepCount++;
    }
    // bw path
    memset(slaveTrajBw, 0., dimension * eSSdim * sizeof(double));
    for (int i = 0; i < eSSdim; ++i)
        slaveTrajBw[i * dimension + i] += eps;

    for (int t = 0; t < Tnorm; t += 1) {
        --stepCount;
        for (int i = 0; i < eSSdim; i++) {
            dverkStepVarMat(&slaveTrajBw[i * dimension], dimension, diffFuncVarRev, params, step, arg,
                traj_dots + stepCount * dimension, k1, k2, k3, k4, k5, k6, k7, k8);
        }
        ortVecs(slaveTrajBw, dimension, eSSdim, projSum);
        normalizeVecs(slaveTrajBw, dimension, eSSdim, eps);
    }
    for (int t = 0; t < T; t += 1) {
        stepCount--;
        for (int i = 0; i < eSSdim; i++)
            dverkStepVarMat(&slaveTrajBw[i * dimension], dimension, diffFuncVarRev, params, step, arg,
                traj_dots + stepCount * dimension, k1, k2, k3, k4, k5, k6, k7, k8);

        ortVecs(slaveTrajBw, dimension, eSSdim, projSum);

        for (int i = 0; i < eSSdim; i++)
            lyapExpBw[i] += log(vecNorm(&slaveTrajBw[i * dimension], dimension) / eps);


        normalizeVecs(slaveTrajBw, dimension, eSSdim, eps);
        for (int i = 0; i < eSSdim; i++) {
            for (int j = 0; j < eSSdim; j++) {
                dotProductsMatrix[i*eSSdim + j] = dotProduct(ncu_s + stepCount * dimension * eSSdim + i * dimension, slaveTrajBw + j * dimension, dimension);
                //printf("%lf ", dotProduct(ncu_s + stepCount * dimension * eSSdim + i * dimension, slaveTrajBw + j * dimension, dimension));
            }
            //std::cout << std::endl;
        }
        //std::cout << std::endl;
        memcpy(dotProductsMatrixT, dotProductsMatrix, eSSdim * eSSdim * sizeof(double));
        matrixTranspose(dotProductsMatrixT, eSSdim, eSSdim);
        matrixMult(A, dotProductsMatrixT, eSSdim, eSSdim, dotProductsMatrix, eSSdim, eSSdim);
        findEigenValues(A, R, Q, projSum, eSSdim, eSSdim, 20, eigenValues);
        compare(eigenValues, eSSdim);
        double minSingularVal = pow(eigenValues[0], 0.5);
        double phi_cur = std::abs(3.14159265358979323846264338327950288 / 2. - acos(minSingularVal));
        if (phi_cur < *phi) {
            *phi = phi_cur;
        }
    }
    for (int i = 0; i < expsNum; i++)
        lyapExpBw[i] = lyapExpBw[i] / calcTime;
    return 0;
}

void cudaAngleNode2Warm(double *main, double *sub, double *essCur, const int dimension,
                       void(*diffFunc)(const double*, double*, const double*), void(*diffFuncVar)(const double*, double*, const double*, const double*),
                       const double *params, double eps, int expsNum, double step, double timeSkip, double timeSkipClP, const int eSSdim) {
    double projSum[MAX_DIMENSION], arg[MAX_DIMENSION];
    double k1[MAX_DIMENSION], k2[MAX_DIMENSION], k3[MAX_DIMENSION], k4[MAX_DIMENSION], k5[MAX_DIMENSION], k6[MAX_DIMENSION], k7[MAX_DIMENSION], k8[MAX_DIMENSION];
    memset(k1, 0., MAX_DIMENSION * sizeof(double));
    memset(k2, 0., MAX_DIMENSION * sizeof(double));
    memset(k3, 0., MAX_DIMENSION * sizeof(double));
    memset(k4, 0., MAX_DIMENSION * sizeof(double));
    memset(k5, 0., MAX_DIMENSION * sizeof(double));
    memset(k6, 0., MAX_DIMENSION * sizeof(double));
    memset(k7, 0., MAX_DIMENSION * sizeof(double));
    memset(k8, 0., MAX_DIMENSION * sizeof(double));
    memset(arg, 0., MAX_DIMENSION * sizeof(double));
    memset(projSum, 0., MAX_DIMENSION * sizeof(double));
    double Twarm = timeSkip / step;
    double Tnorm = timeSkipClP / step;
    //printf("TIME FOR SKIP: %lf\n", Twarm);
    //printf("TIME FOR NORW: %lf\n", Tnorm);
    for (int t = 0; t < Twarm; t += 1) // Skip points to reach the attractor
        dverkStep(main, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);

    memset(sub, 0., expsNum * dimension * sizeof(double));

    if (expsNum < dimension - eSSdim)
        expsNum = dimension - eSSdim;
    for (int i = 0; i < expsNum; i++) {
        sub[i * dimension + i] += eps;
    }
    for (int t = 0; t < Tnorm; t += 1) {
        for (int i = 0; i < expsNum; i++)
            dverkStepVarMat(&sub[i * dimension], dimension, diffFuncVar,
                params, step, arg, main, k1, k2, k3, k4, k5, k6, k7, k8);
        dverkStep(main, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        ortVecs(sub, dimension, expsNum, projSum);
        normalizeVecs(sub, dimension, expsNum, eps);
    }

}

int cudaAngleNode2Calculate(double *mainTrajectory, double *slaveTrajectories, double *essCur, const int dimension,
                        void(*diffFunc)(const double*, double*,const double*), void(*diffFuncVar)(const double*, double*,const double*,const double*),
                        void(*diffFuncVarTrans)(const double*, double*,const double*,const double*),
                        const double *params, double eps, double *lyapExp, double *lyapExpBw, double *phi,
                        int expsNum, double step, double timeSkip, double timeSkipClP, double calcTime,
                        double addSkipSteps, double *traj_dots, double *ncu_s, const int eSSdim) {
    double projSum[MAX_DIMENSION], arg[MAX_DIMENSION];

    double main_afterForward[MAX_DIMENSION]; // to calc main traj after fw time is completed
    double slaveTrajBw[MAX_DIMENSION * MAX_DIMENSION]; // exp_num * max_dimension
    double dotProductsMatrix[MAX_DIMENSION * MAX_DIMENSION];
    double dotProductsMatrixT[MAX_DIMENSION * MAX_DIMENSION];
    double A[MAX_DIMENSION * MAX_DIMENSION];
    double R[MAX_DIMENSION * MAX_DIMENSION];
    double Q[MAX_DIMENSION * MAX_DIMENSION];
    double eigenValues[MAX_DIMENSION];
    memset(dotProductsMatrix, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(dotProductsMatrixT, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(A, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(R, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(Q, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(eigenValues, 0., MAX_DIMENSION * sizeof(double));

    memset(main_afterForward, 0., MAX_DIMENSION * sizeof(double));
    memset(lyapExp, 0., MAX_DIMENSION * sizeof(double));
    memset(lyapExpBw, 0., MAX_DIMENSION * sizeof(double));

    if (expsNum < dimension - eSSdim)
        expsNum = dimension - eSSdim;
    int eCUdim = dimension - eSSdim;
    //std::cout << "CU: " << eCUdim << std::endl;

    int stepCount = 0;
    double k1[MAX_DIMENSION], k2[MAX_DIMENSION], k3[MAX_DIMENSION], k4[MAX_DIMENSION], k5[MAX_DIMENSION], k6[MAX_DIMENSION], k7[MAX_DIMENSION], k8[MAX_DIMENSION];
    memset(k1, 0., MAX_DIMENSION * sizeof(double));
    memset(k2, 0., MAX_DIMENSION * sizeof(double));
    memset(k3, 0., MAX_DIMENSION * sizeof(double));
    memset(k4, 0., MAX_DIMENSION * sizeof(double));
    memset(k5, 0., MAX_DIMENSION * sizeof(double));
    memset(k6, 0., MAX_DIMENSION * sizeof(double));
    memset(k7, 0., MAX_DIMENSION * sizeof(double));
    memset(k8, 0., MAX_DIMENSION * sizeof(double));
    memset(arg, 0., MAX_DIMENSION * sizeof(double));
    memset(projSum, 0., MAX_DIMENSION * sizeof(double));
    double T = calcTime / step;
    double Tnorm = timeSkipClP / step;
    for (int t = 0; t < T; t += 1) {
        for (int i = 0; i < expsNum; i++)
            dverkStepVarMat(&slaveTrajectories[i * dimension], dimension, diffFuncVar, params, step, arg,
                mainTrajectory, k1, k2, k3, k4, k5, k6, k7, k8);

        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        
        ortVecs(slaveTrajectories, dimension, expsNum, projSum);

        for (int i = 0; i < expsNum; i++)
            lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / eps);

        normalizeVecs(slaveTrajectories, dimension, expsNum, eps);
        memcpy(ncu_s + stepCount * dimension * eCUdim, slaveTrajectories, sizeof(double) * dimension * eCUdim);
        memcpy(traj_dots + stepCount * dimension, mainTrajectory, sizeof(double) * dimension);
        stepCount++;
    }

    for (int i = 0; i < expsNum; i++) {
        lyapExp[i] = lyapExp[i] / calcTime;
    }

    memcpy(main_afterForward, mainTrajectory, sizeof(double) * dimension);
    for (int t = 0; t < Tnorm; t += 1) {
        dverkStep(main_afterForward, dimension, diffFunc,
            params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        memcpy(traj_dots + stepCount * dimension, main_afterForward, sizeof(double) * dimension);
        stepCount++;
    }
    // bw path
    memset(slaveTrajBw, 0., dimension * eCUdim * sizeof(double));
    for (int i = 0; i < eCUdim; ++i)
        slaveTrajBw[i * dimension + i] += eps;

    for (int t = 0; t < Tnorm; t += 1) {
        --stepCount;
        for (int i = 0; i < eCUdim; i++) {
            dverkStepVarMat(&slaveTrajBw[i * dimension], dimension, diffFuncVarTrans, params, step, arg,
                traj_dots + stepCount * dimension, k1, k2, k3, k4, k5, k6, k7, k8);
        }
        ortVecs(slaveTrajBw, dimension, eCUdim, projSum);
        normalizeVecs(slaveTrajBw, dimension, eCUdim, eps);
    }

    for (int t = 0; t < T; t += 1) {
        stepCount--;
        for (int i = 0; i < eCUdim; i++)
            dverkStepVarMat(&slaveTrajBw[i * dimension], dimension, diffFuncVarTrans, params, step, arg,
                traj_dots + stepCount * dimension, k1, k2, k3, k4, k5, k6, k7, k8);

        ortVecs(slaveTrajBw, dimension, eCUdim, projSum);

        for (int i = 0; i < eCUdim; i++)
            lyapExpBw[i] += log(vecNorm(&slaveTrajBw[i * dimension], dimension) / eps);


        normalizeVecs(slaveTrajBw, dimension, eCUdim, eps);
        for (int i = 0; i < eCUdim; i++)
            for (int j = 0; j < eCUdim; j++) {
                dotProductsMatrix[i*eCUdim + j] = dotProduct(ncu_s + stepCount * dimension * eCUdim + i * dimension, slaveTrajBw + j * dimension, dimension);
            }
        memcpy(dotProductsMatrixT, dotProductsMatrix, eCUdim * eCUdim * sizeof(double));
        matrixTranspose(dotProductsMatrixT, eCUdim, eCUdim);
        matrixMult(A, dotProductsMatrixT, eCUdim, eCUdim, dotProductsMatrix, eCUdim, eCUdim);
        findEigenValues(A, R, Q, projSum, eCUdim, eCUdim, 50, eigenValues);
        compare(eigenValues, eCUdim);
        double minSingularVal = pow(eigenValues[0], 0.5);
        double phi_cur = std::abs(3.14159265358979323846264338327950288 / 2. - acos(minSingularVal));
        if (phi_cur < *phi) {
            *phi = phi_cur;
        }
    }
    for (int i = 0; i < eCUdim; i++)
        lyapExpBw[i] = lyapExpBw[i] / calcTime;
    return 0;
}