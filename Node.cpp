#include "consts.h"
#include "integrator.cpp"
#include "iterator.cpp"
#include "linal.cpp"
#include "utils.h"
#include <cstring>
#include <stdlib.h>
#include <cstdlib>


void cudaAngleCUNodeWarm(double *main, double *sub, double *essCur, const int dimension,
                         void(*diffFunc)(const double*, double*, const double*), void(*diffFuncVar)(const double*, double*, const double*, const double*),
                         double *params, double eps, int expsNum, double step, double timeSkip, 
                         double timeSkipClP, double addSkipSteps, const bool type, int essDim, MemoryForIntegrator * memI) {
     /*Выход траектории на аттрактор и прогрев векторов вдоль траектории для линеаризованной системы*/
    //Будем считать столько показателей Ляпунова сколько необходимо для алгоритма или больше, если задал пользователь
    if (expsNum<essDim)
        expsNum = essDim;

    int skipCount = 0;
    if (type) {
        int Twarm = timeSkip / step;
        int Tnorm = timeSkipClP / step;

        for (int t = 0; t < Twarm; t += 1) // Skip points to reach the attractor
            dverkStep(main, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

        memset(sub, 0., expsNum * dimension * sizeof(double));

        for (int i = 0; i < expsNum; i++)
            sub[i * dimension + i] += 1;

        for (int t = 0; t < Tnorm; t += 1) {
            for (int i = 0; i < expsNum; i++)
                dverkStepVarMat(&sub[i * dimension], dimension, diffFuncVar,
                                params, step, memI->arg, main, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            dverkStep(main, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            if (addSkipSteps == skipCount) {
                ortVecs(sub, dimension, expsNum, memI->projSum);
                normalizeVecs(sub, dimension, expsNum, 1);
                skipCount = 0;
            }
            skipCount++;
        }
    }
    else {
        for (double t = 0; t < timeSkip; t += 1) // Skip points to reach the attractor
            StepMAP(main, diffFunc, params, memI->arg, dimension);

        memset(sub, 0., expsNum * dimension * sizeof(double));

        for (int i = 0; i < expsNum; i++)
            sub[i * dimension + i] += eps;

        for (int t = 0; t < timeSkipClP; t += 1) {
            for (int i = 0; i < expsNum; i++)
                StepMAPVAR(&sub[i * dimension], diffFuncVar, params, memI->arg, dimension, main);
            StepMAP(main, diffFunc, params, memI->arg, dimension);
            if (addSkipSteps == skipCount) {
                ortVecs(sub, dimension, expsNum, memI->projSum);
                normalizeVecs(sub, dimension, expsNum, 1);
                skipCount = 0;
            }
            skipCount++;
        }
    }

}


int cudaAngleCUNodeCalculate(double *mainTrajectory, double *slaveTrajectories, double *essCur, const int dimension,
                             void(*diffFunc)(const double*, double*, const double*), void(*diffFuncVar)(const double*, double*, const double*, const double*), 
                             void(*diffFuncVarTrans)(const double*, double*,const double*,const double*),
                             double *params, double eps, double *lyapExp, double *lyapExpBw, double *phi,
                             int expsNum, double step, double timeSkip, double timeSkipClP, double calcTime,
                             double addSkipSteps, double *traj_dots, double *ncu_s, const bool type, int essDim, MemoryForIntegrator * memI, MemoryForMatrix * memM, int eigenIteration) {
    /*Алгоритм в случае когда размерность E^{cu} меньше или равна половине от размерности системы*/
    if (expsNum < essDim)
        expsNum = essDim;

    double main_afterForward[MAX_DIMENSION]; // to calc main traj after fw time is completed

    memset(main_afterForward, 0., MAX_DIMENSION * sizeof(double));
    memset(lyapExp, 0., expsNum * sizeof(double));
    memset(lyapExpBw, 0., expsNum * sizeof(double));

    size_t T = calcTime / step;
    size_t Tnorm = timeSkipClP / step;

    int skipCount = 0;
    int stepCountPoints = 0;
    int stepCountVectors = 0;
    if (type) {
        for (int t = 0; t < T; t += 1) {
            skipCount++;
            for (int i = 0; i < expsNum; i++) {
                dverkStepVarMat(&slaveTrajectories[i * dimension], dimension, diffFuncVar, params, step, memI->arg,
                                mainTrajectory, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            }

            dverkStep(mainTrajectory, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajectories, dimension, expsNum, memI->projSum);

                for (int i = 0; i < expsNum; i++)
                    lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / 1);

                normalizeVecs(slaveTrajectories, dimension, expsNum, 1);
                memcpy(ncu_s + stepCountVectors * dimension * essDim, slaveTrajectories, sizeof(double) * dimension * essDim);
                stepCountVectors++;
                skipCount = 0;
            }
            memcpy(traj_dots + stepCountPoints * dimension, mainTrajectory, sizeof(double) * dimension);
            stepCountPoints++;
        }

        for (int i = 0; i < expsNum; i++) {
            lyapExp[i] = lyapExp[i] / calcTime;
        }

        memcpy(main_afterForward, mainTrajectory, sizeof(double) * dimension);
        for (int t = 0; t < Tnorm; t += 1) {
            dverkStep(main_afterForward, dimension, diffFunc,
                      params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            memcpy(traj_dots + stepCountPoints * dimension, main_afterForward, sizeof(double) * dimension);
            stepCountPoints++;
        }
        // bw path
        double slaveTrajBw[MAX_DIMENSION * MAX_DIMENSION]; // exp_num * max_dimension
        memset(slaveTrajBw, 0., dimension * essDim * sizeof(double));
        for (int i = 0; i < essDim; ++i) {
            slaveTrajBw[i * dimension + i] += 1;
        }
        skipCount = 0;
        for (int t = 0; t < Tnorm; t += 1) {
            --stepCountPoints;
            for (int i = 0; i < essDim; i++)
                dverkStepVarMat(&slaveTrajBw[i*dimension], dimension, diffFuncVarTrans, params, step, memI->arg,
                                traj_dots + stepCountPoints * dimension, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                skipCount = 0;
            }
            skipCount++;
        }

        skipCount = addSkipSteps;
        for (int t = 0; t < T; t += 1) {
            stepCountPoints--;
            for (int i = 0; i < essDim; i++)
                dverkStepVarMat(&slaveTrajBw[i*dimension], dimension, diffFuncVarTrans, params, step, memI->arg,
                                traj_dots + stepCountPoints * dimension, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            if (skipCount == addSkipSteps) {
                stepCountVectors--;
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                for (int j = 0; j < essDim; j++)
                    for (int i = 0; i < essDim; i++) {
                        memM->dotProductsMatrix[j*essDim + i] = dotProduct(ncu_s + stepCountVectors * dimension * essDim + j * dimension, slaveTrajBw + i * dimension, dimension);
                    }
                memcpy(memM->dotProductsMatrixT, memM->dotProductsMatrix, essDim * essDim * sizeof(double));
                matrixTranspose(memM->dotProductsMatrixT, essDim, essDim);
                matrixMult(memM->A, memM->dotProductsMatrixT, essDim, essDim, memM->dotProductsMatrix, essDim, essDim);
                double angle = pow(findEigenValues(memM->A, memM->R, memM->Q, memI->projSum, essDim, essDim, 50, memM->eigenValues), 0.5);
                double phi_cur;
                phi_cur = std::abs(3.14159265358979323846264338327950288 / 2. - acos(angle));
                if (phi_cur < *phi) {
                    *phi = phi_cur;
                }
                skipCount = 0;
            }
            skipCount++;
        }
    }
    else {
        for (int t = 0; t < calcTime; t += 1) {
            for (int i = 0; i < expsNum; i++)
                StepMAPVAR(&slaveTrajectories[i * dimension], diffFuncVar, params, memI->arg, dimension, mainTrajectory);

            StepMAP(mainTrajectory, diffFunc, params, memI->arg, dimension);

            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajectories, dimension, expsNum, memI->projSum);

                for (int i = 0; i < expsNum; i++)
                    lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / 1);

                normalizeVecs(slaveTrajectories, dimension, expsNum, 1);
                memcpy(ncu_s + stepCountVectors * dimension * essDim, slaveTrajectories, sizeof(double) * dimension * essDim);
                stepCountVectors++;
                skipCount = 0;
            }
            memcpy(traj_dots + stepCountPoints * dimension, mainTrajectory, sizeof(double) * dimension);
            stepCountPoints++;
        }

        for (int i = 0; i < expsNum; i++) {
            lyapExp[i] = lyapExp[i] / calcTime;
        }

        memcpy(main_afterForward, mainTrajectory, sizeof(double) * dimension);
        for (int t = 0; t < timeSkipClP; t += 1) {
            StepMAP(main_afterForward, diffFunc, params, memI->arg, dimension);
            memcpy(traj_dots + stepCountPoints * dimension, main_afterForward, sizeof(double) * dimension);
            stepCountPoints++;
        }
        // bw path
        double slaveTrajBw[MAX_DIMENSION * MAX_DIMENSION]; // exp_num * max_dimension
        memset(slaveTrajBw, 0., dimension * essDim * sizeof(double));
        for (int i = 0; i < essDim; ++i) {
            slaveTrajBw[i * dimension + i] += 1;
        }

        skipCount = 0;
        for (int t = 0; t < timeSkipClP; t += 1) {
            --stepCountPoints;
            for (int i = 0; i < essDim; i++)
                StepMAPVAR(&slaveTrajBw[i*dimension], diffFuncVarTrans, params, memI->arg, dimension, traj_dots + stepCountPoints * dimension);
            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, expsNum, 1);
                skipCount = 0;
            }
            skipCount++;
        }
        for (int t = 0; t < calcTime; t += 1) {
            stepCountPoints--;
            for (int i = 0; i < essDim; i++)
                StepMAPVAR(&slaveTrajBw[i*dimension], diffFuncVarTrans, params, memI->arg, dimension, traj_dots + stepCountPoints * dimension);
            if (skipCount == addSkipSteps) {
                stepCountVectors--;
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                for (int j = 0; j < essDim; j++)
                    for (int i = 0; i < essDim; i++) {
                        memM->dotProductsMatrix[j*essDim + i] = dotProduct(ncu_s + stepCountVectors * dimension * essDim + j * dimension, slaveTrajBw + i * dimension, dimension);
                    }
                memcpy(memM->dotProductsMatrixT, memM->dotProductsMatrix, essDim * essDim * sizeof(double));
                matrixTranspose(memM->dotProductsMatrixT, essDim, essDim);
                matrixMult(memM->A, memM->dotProductsMatrixT, essDim, essDim, memM->dotProductsMatrix, essDim, essDim);
                double angle = pow(findEigenValues(memM->A, memM->R, memM->Q, memI->projSum, essDim, essDim, 50, memM->eigenValues), 0.5);
                double phi_cur;
                phi_cur = std::abs(3.14159265358979323846264338327950288 / 2. - acos(angle));
                if (phi_cur < *phi) {
                    *phi = phi_cur;
                }
                skipCount = 0;
            }
            skipCount++;
        }
    }
    return 0;
}


void cudaAngleSSNodeWarm(double *main, double *sub, double *essCur, const int dimension,
                         void(*diffFunc)(const double*, double*, const double*), void(*diffFuncVar)(const double*, double*, const double*, const double*), void(*diffFuncVarTransRev)(const double*, double*,const double*,const double*),
                         double *params, double eps, int expsNum, double step, double timeSkip, double timeSkipClP, const bool type, const int essDim, MemoryForIntegrator * memI, double addSkipSteps) {
    /*Выход траектории на аттрактор и прогрев векторов вдоль траектории для транспонированной линеаризованной обратной системы*/
    
    int skipCount = 0;
    if (type) {
        int Twarm = timeSkip / step;
        int Tnorm = timeSkipClP / step;

        for (int t = 0; t < Twarm; t += 1) // Skip points to reach the attractor
            dverkStep(main, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

        memset(sub, 0., expsNum * dimension * sizeof(double));
        memset(essCur, 0., essDim * dimension * sizeof(double));

        for (int i = 0; i < expsNum; i++)
            sub[i * dimension + i] += 1;

        for (int i = 0; i < essDim; i++)
            essCur[i * dimension + i] += 1;

        for (int t = 0; t < Tnorm; t += 1) {
            skipCount++;
            for (int i = 0; i < expsNum; i++)
                dverkStepVarMat(&sub[i * dimension], dimension, diffFuncVar,
                    params, step, memI->arg, main, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

            for (int j = 0; j < essDim; j++)
                dverkStepVarMat(&essCur[j*dimension], dimension, diffFuncVarTransRev, params, step, memI->arg, main, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

            dverkStep(main, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            if (addSkipSteps == skipCount) {
                ortVecs(sub, dimension, expsNum, memI->projSum);
                ortVecs(essCur, dimension, essDim, memI->projSum);
                normalizeVecs(sub, dimension, expsNum, 1);
                normalizeVecs(essCur, dimension, essDim, 1);
                skipCount = 0;
            }
        }
    }
    else {
        for (double t = 0; t < timeSkip; t += 1) // Skip points to reach the attractor
            StepMAP(main, diffFunc, params, memI->arg, dimension);

        memset(sub, 0., expsNum * dimension * sizeof(double));
        memset(essCur, 0., essDim * dimension * sizeof(double));

        for (int i = 0; i < expsNum; i++)
            sub[i * dimension + i] += 1;
        
        for (int i = 0; i<essDim; i++)
            essCur[i * dimension + i] += 1;
        

        for (int t = 0; t < timeSkipClP; t += 1) {
            skipCount++;
            for (int i = 0; i < expsNum; i++)
                StepMAPVAR(&sub[i * dimension], diffFuncVar, params, memI->arg, dimension, main);
            for (int j = 0; j < essDim; j++)
                StepMAPVAR(&essCur[j*dimension], diffFuncVarTransRev, params, memI->arg, dimension, main);

            StepMAP(main, diffFunc, params, memI->arg, dimension);
            if (addSkipSteps == skipCount) {
                ortVecs(sub, dimension, expsNum, memI->projSum);
                ortVecs(essCur, dimension, essDim, memI->projSum);
                normalizeVecs(sub, dimension, expsNum, 1);
                normalizeVecs(essCur, dimension, essDim, 1);
                skipCount = 0;
            }
        }
    }

}


int cudaAngleSSNodeCalculate(double *mainTrajectory, double *slaveTrajectories, double *essCur, const int dimension,
                             void(*diffFunc)(const double*, double*, const double*), void(*diffFuncVar)(const double*, double*, const double*, const double*),
                             void(*diffFuncVarRev)(const double*, double*,const double*,const double*), 
                             void(*diffFuncVarTransRev)(const double*, double*,const double*,const double*),
                             double *params, double eps, double *lyapExp, double *lyapExpBw, double *phi,
                             int expsNum, double step, double timeSkip, double timeSkipClP, double calcTime,
                             double addSkipSteps, double *traj_dots, double *ncu_s, const bool type, const int essDim,
                             MemoryForIntegrator * memI, MemoryForMatrix * memM, int eigenIteration) {

    /*Алгоритм в случае когда размерность E^{cu} больше половины от размерности системы*/
    double main_afterForward[MAX_DIMENSION]; // to calc main traj after fw time is completed

    memset(main_afterForward, 0., MAX_DIMENSION * sizeof(double));
    memset(lyapExp, 0., expsNum * sizeof(double));
    memset(lyapExpBw, 0., expsNum * sizeof(double));

    size_t T = calcTime / step;
    size_t Tnorm = timeSkipClP / step;

    int skipCount = 0;
    int stepCountPoints = 0;
    int stepCountVectors = 0;
    if (type) {
        for (int t = 0; t < T; t += 1) {
            skipCount++;
            for (int i = 0; i < expsNum; i++)
                dverkStepVarMat(&slaveTrajectories[i * dimension], dimension, diffFuncVar, params, step, memI->arg,
                    mainTrajectory, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

            for (int i = 0; i < essDim; i++)
                dverkStepVarMat(&essCur[i * dimension], dimension, diffFuncVarTransRev, params, step, memI->arg,
                                mainTrajectory, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

            dverkStep(mainTrajectory, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajectories, dimension, expsNum, memI->projSum);
                ortVecs(essCur, dimension, essDim, memI->projSum);

                for (int i = 0; i < expsNum; i++)
                    lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / 1);

                normalizeVecs(slaveTrajectories, dimension, expsNum, 1);
                normalizeVecs(essCur, dimension, essDim, 1);
                memcpy(ncu_s + stepCountVectors * dimension * essDim, essCur, sizeof(double) * dimension * essDim);
                stepCountVectors++;
                skipCount = 0;
            }
            memcpy(traj_dots + stepCountPoints * dimension, mainTrajectory, sizeof(double) * dimension);
            stepCountPoints++;
        }

        for (int i = 0; i < expsNum; i++) {
            lyapExp[i] = lyapExp[i] / calcTime;
        }

        memcpy(main_afterForward, mainTrajectory, sizeof(double) * dimension);
        for (int t = 0; t < Tnorm; t += 1) {
            dverkStep(main_afterForward, dimension, diffFunc,
                params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            memcpy(traj_dots + stepCountPoints * dimension, main_afterForward, sizeof(double) * dimension);
            stepCountPoints++;
        }
        // bw path
        double slaveTrajBw[MAX_DIMENSION * MAX_DIMENSION]; // exp_num * max_dimension
        memset(slaveTrajBw, 0., dimension * essDim * sizeof(double));
        for (int i = 0; i < essDim; ++i) {
            slaveTrajBw[i * dimension + i] += 1;
        }

        skipCount = 0;
        for (int t = 0; t < Tnorm; t += 1) {
            --stepCountPoints;
            for (int i = 0; i < essDim; i++)
                dverkStepVarMat(&slaveTrajBw[i*dimension], dimension, diffFuncVarRev, params, step, memI->arg,
                    traj_dots + stepCountPoints * dimension, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                skipCount=0;
            }
            skipCount++;
        }

        skipCount = addSkipSteps;
        for (int t = 0; t < T; t += 1) {
            stepCountPoints--;
            for (int i = 0; i < essDim; i++)
                dverkStepVarMat(&slaveTrajBw[i*dimension], dimension, diffFuncVarRev, params, step, memI->arg,
                    traj_dots + stepCountPoints * dimension, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            if (skipCount == addSkipSteps) {
                stepCountVectors--;
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                for (int j = 0; j < essDim; j++)
                    for (int i = 0; i < essDim; i++)
                        memM->dotProductsMatrix[j*essDim + i] = dotProduct(ncu_s + stepCountVectors * dimension * essDim + j * dimension, slaveTrajBw + i * dimension, dimension);
                memcpy(memM->dotProductsMatrixT, memM->dotProductsMatrix, essDim * essDim * sizeof(double));
                matrixTranspose(memM->dotProductsMatrixT, essDim, essDim);
                matrixMult(memM->A, memM->dotProductsMatrixT, essDim, essDim, memM->dotProductsMatrix, essDim, essDim);
                double angle = pow(findEigenValues(memM->A, memM->R, memM->Q, memI->projSum, essDim, essDim, eigenIteration, memM->eigenValues), 0.5);
                double phi_cur;
                phi_cur = std::abs(3.14159265358979323846264338327950288 / 2. - acos(angle));
                if (phi_cur < *phi) {
                    *phi = phi_cur;
                }
                skipCount = 0;
            }
            skipCount++;
        }
    }
    else {
        for (int t = 0; t < calcTime; t += 1) {
            skipCount++;
            for (int i = 0; i < expsNum; i++)
                StepMAPVAR(&slaveTrajectories[i * dimension], diffFuncVar, params, memI->arg, dimension, mainTrajectory);
            for (int i = 0; i < essDim; i++)
                StepMAPVAR(&essCur[i*dimension], diffFuncVarTransRev, params, memI->arg, dimension, mainTrajectory);

            StepMAP(mainTrajectory, diffFunc, params, memI->arg, dimension);

            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajectories, dimension, expsNum, memI->projSum);
                ortVecs(essCur, dimension, essDim, memI->projSum);

                for (int i = 0; i < expsNum; i++)
                    lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / 1);

                normalizeVecs(slaveTrajectories, dimension, expsNum, 1);
                normalizeVecs(essCur, dimension, essDim, 1);
                memcpy(ncu_s + stepCountVectors * dimension * essDim, essCur, sizeof(double) * dimension * essDim);
                stepCountVectors++;
                skipCount = 0;
            }
            memcpy(traj_dots + stepCountPoints * dimension, mainTrajectory, sizeof(double) * dimension);
            stepCountPoints++;
        }

        for (int i = 0; i < expsNum; i++) {
            lyapExp[i] = lyapExp[i] / calcTime;
        }

        memcpy(main_afterForward, mainTrajectory, sizeof(double) * dimension);
        for (int t = 0; t < timeSkipClP; t += 1) {
            StepMAP(main_afterForward, diffFunc, params, memI->arg, dimension);
            memcpy(traj_dots + stepCountPoints * dimension, main_afterForward, sizeof(double) * dimension);
            stepCountPoints++;
        }
        // bw path
        double slaveTrajBw[MAX_DIMENSION * MAX_DIMENSION]; // exp_num * max_dimension
        memset(slaveTrajBw, 0., dimension * essDim * sizeof(double));
        for (int i = 0; i < essDim; ++i) {
            slaveTrajBw[i * dimension + i] += 1;
        }

        skipCount = 0;
        for (int t = 0; t < timeSkipClP; t += 1) {
            --stepCountPoints;
            for (int i = 0; i < essDim; i++)
                StepMAPVAR(&slaveTrajBw[i*dimension], diffFuncVarRev, params, memI->arg, dimension, traj_dots + stepCountPoints * dimension);
            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                skipCount = 0;
            }
            skipCount++;
        }
        for (int t = 0; t < calcTime; t += 1) {
            stepCountPoints--;
            for (int i = 0; i < essDim; i++)
                StepMAPVAR(&slaveTrajBw[i*dimension], diffFuncVarRev, params, memI->arg, dimension, traj_dots + stepCountPoints * dimension);
            if (skipCount == addSkipSteps) {
                stepCountVectors--;
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                for (int j = 0; j < essDim; j++)
                    for (int i = 0; i < essDim; i++)
                        memM->dotProductsMatrix[j*essDim + i] = dotProduct(ncu_s + stepCountVectors * dimension * essDim + j * dimension, slaveTrajBw + i * dimension, dimension);
                memcpy(memM->dotProductsMatrixT, memM->dotProductsMatrix, essDim * essDim * sizeof(double));
                matrixTranspose(memM->dotProductsMatrixT, essDim, essDim);
                matrixMult(memM->A, memM->dotProductsMatrixT, essDim, essDim, memM->dotProductsMatrix, essDim, essDim);
                double angle = pow(findEigenValues(memM->A, memM->R, memM->Q, memI->projSum, essDim, essDim, eigenIteration, memM->eigenValues), 0.5);
                double phi_cur;
                phi_cur = std::abs(3.14159265358979323846264338327950288 / 2. - acos(angle));
                if (phi_cur < *phi) {
                    *phi = phi_cur;
                }
                skipCount = 0;
            }
            skipCount++;
        }
    }
    return 0;
}

void cudaAngleSlaveNodeWarm(double *main, double *sub, double *essCur, const int dimension,
                         void(*diffFunc)(const double*, double*, const double*),
                         double *params, double eps, int expsNum, double step, double timeSkip,
                         double timeSkipClP, const bool type, const int essDim, MemoryForIntegrator * memI,
                         double addSkipSteps) {
    /*Выход траектории на аттрактор и прогрев векторов вдоль траектории для транспонированной линеаризованной обратной системы*/
    
    int skipCount = 0;
    if (type) {
        int Twarm = timeSkip / step;
        int Tnorm = timeSkipClP / step;

        for (int t = 0; t < Twarm; t += 1) // Skip points to reach the attractor
            dverkStep(main, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);


        for (int i = 0; i < dimension; i++){
            memcpy(&sub[i*dimension], main, dimension * sizeof(double));   
            sub[i * dimension + i] += eps;
        }

        for (int t = 0; t < Tnorm; t += 1) {
            skipCount++;
            dverkStep(main, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            
            for (int i = 0; i < dimension; i++) {
                dverkStep(&sub[i * dimension], dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
                for (int j = 0; j < dimension; j++)
                    sub[i * dimension + j] -= main[j];
            }
            
            if (addSkipSteps == skipCount) {
                ortVecs(sub, dimension, dimension, memI->projSum);
                normalizeVecs(sub, dimension, dimension, eps);
                skipCount = 0;
            }

            for (int i = 0; i < expsNum; i++) {
                for (int j = 0; j < dimension; j++)
                    sub[i * dimension + j] += main[j];
            }
        }
    }
    else {
        for (double t = 0; t < timeSkip; t += 1) // Skip points to reach the attractor
            StepMAP(main, diffFunc, params, memI->arg, dimension);

        memset(sub, 0., expsNum * dimension * sizeof(double));
        memset(essCur, 0., essDim * dimension * sizeof(double));

        for (int i = 0; i < expsNum; i++)
            sub[i * dimension + i] += 1;
        
        for (int i = 0; i<essDim; i++)
            essCur[i * dimension + i] += 1;
        

        for (int t = 0; t < timeSkipClP; t += 1) {
            skipCount++;
            StepMAP(main, diffFunc, params, memI->arg, dimension);
            if (addSkipSteps == skipCount) {
                ortVecs(sub, dimension, expsNum, memI->projSum);
                ortVecs(essCur, dimension, essDim, memI->projSum);
                normalizeVecs(sub, dimension, expsNum, 1);
                normalizeVecs(essCur, dimension, essDim, 1);
                skipCount = 0;
            }
        }
    }

}


int cudaAngleSlaveNodeCalculate(double *mainTrajectory, double *slaveTrajectories, double *essCur, const int dimension,
                             void(*diffFunc)(const double*, double*, const double*),
                             void(*diffFuncRev)(const double*, double*,const double*), 
                             double *params, double eps, double *lyapExp, double *lyapExpBw, double *phi,
                             int expsNum, double step, double timeSkip, double timeSkipClP, double calcTime,
                             double addSkipSteps, double *traj_dots, double *ncu_s, const bool type, const int essDim,
                             MemoryForIntegrator * memI, MemoryForMatrix * memM, int eigenIteration) {
    /*Алгоритм в случае когда размерность E^{cu} больше половины от размерности системы*/
    double main_afterForward[MAX_DIMENSION]; // to calc main traj after fw time is completed

    memset(main_afterForward, 0., MAX_DIMENSION * sizeof(double));
    memset(lyapExp, 0., expsNum * sizeof(double));
    memset(lyapExpBw, 0., expsNum * sizeof(double));

    size_t T = calcTime / step;
    size_t Tnorm = timeSkipClP / step;

    int skipCount = 0;
    int stepCountPoints = 0;
    int stepCountVectors = 0;
    if (type) {
        for (int t = 0; t < T; t += 1) {
            skipCount++;
            dverkStep(mainTrajectory, dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);

            for (int i = 0; i < dimension; i++) {
                dverkStep(&slaveTrajectories[i * dimension], dimension, diffFunc, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
                for (int j = 0; j < dimension; j++)
                    slaveTrajectories[i * dimension + j] -= mainTrajectory[j];
            }

            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajectories, dimension, dimension, memI->projSum);

                for (int i = 0; i < expsNum; i++)
                    lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / eps);

                normalizeVecs(slaveTrajectories, dimension, dimension, 1);
                memcpy(ncu_s + stepCountVectors * dimension * essDim, slaveTrajectories + dimension*(dimension-essDim), sizeof(double) * dimension * essDim);
                
                normalizeVecs(slaveTrajectories, dimension, dimension, eps);
                stepCountVectors++;
                skipCount = 0;
            }
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++)
                    slaveTrajectories[i * dimension + j] += mainTrajectory[j];
            }
            memcpy(traj_dots + stepCountPoints * dimension, mainTrajectory, sizeof(double) * dimension);
            stepCountPoints++;
        }

        for (int i = 0; i < dimension; i++) {
            lyapExp[i] = lyapExp[i] / calcTime;
        }

        memcpy(main_afterForward, mainTrajectory, sizeof(double) * dimension);
        for (int t = 0; t < Tnorm; t += 1) {
            dverkStep(main_afterForward, dimension, diffFunc,
                params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            memcpy(traj_dots + stepCountPoints * dimension, main_afterForward, sizeof(double) * dimension);
            stepCountPoints++;
        }
        // bw path
        double slaveTrajBw[MAX_DIMENSION * MAX_DIMENSION]; // exp_num * max_dimension
        for (int i = 0; i < essDim; ++i) {
            memcpy(&slaveTrajBw[i * dimension], traj_dots + stepCountPoints*dimension, dimension * sizeof(double));
            slaveTrajBw[i * dimension + i] += eps;
        }

        skipCount = 0;
        for (int t = 0; t < Tnorm; t += 1) {
            --stepCountPoints;
            for (int i = 0; i < essDim; i++) {
                dverkStep(&slaveTrajBw[i * dimension], dimension, diffFuncRev, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
                for (int j = 0; j < dimension; j++)
                    slaveTrajBw[i * dimension + j] -= (traj_dots+stepCountPoints*dimension - dimension)[j];
            }
            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, eps);
                skipCount=0;
            }
            for (int i = 0; i < essDim; i++) {
                for (int j = 0; j < dimension; j++)
                    slaveTrajBw[i * dimension + j] += (traj_dots+stepCountPoints*dimension-dimension)[j];
            }
            skipCount++;
        }

        skipCount = addSkipSteps;
        for (int t = 0; t < T; t += 1) {
            stepCountPoints--;
            for (int i = 0; i < essDim; i++) {
                dverkStep(&slaveTrajBw[i * dimension], dimension, diffFuncRev, params, step, memI->arg, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
                for (int j = 0; j < dimension; j++)
                    slaveTrajBw[i * dimension + j] -= (traj_dots+stepCountPoints*dimension - dimension)[j];
            }
            // for (int i = 0; i < essDim; i++)
            //     dverkStepVarMat(&slaveTrajBw[i*dimension], dimension, diffFuncVarRev, params, step, memI->arg,
            //         traj_dots + stepCountPoints * dimension, memI->k1, memI->k2, memI->k3, memI->k4, memI->k5, memI->k6, memI->k7, memI->k8);
            if (skipCount == addSkipSteps) {
                stepCountVectors--;
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                for (int i = 0; i < essDim; i++)
                    lyapExpBw[i] += log(vecNorm(&slaveTrajBw[i * dimension], dimension) / eps);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                for (int j = 0; j < essDim; j++)
                    for (int i = 0; i < essDim; i++)
                        memM->dotProductsMatrix[j*essDim + i] = dotProduct(ncu_s + stepCountVectors * dimension * essDim + j * dimension, slaveTrajBw + i * dimension, dimension);
                memcpy(memM->dotProductsMatrixT, memM->dotProductsMatrix, essDim * essDim * sizeof(double));
                matrixTranspose(memM->dotProductsMatrixT, essDim, essDim);
                matrixMult(memM->A, memM->dotProductsMatrixT, essDim, essDim, memM->dotProductsMatrix, essDim, essDim);
                double angle = pow(findEigenValues(memM->A, memM->R, memM->Q, memI->projSum, essDim, essDim, eigenIteration, memM->eigenValues), 0.5);
                double phi_cur;
                phi_cur = std::abs(3.14159265358979323846264338327950288 / 2. - acos(angle));
                if (phi_cur < *phi) {
                    *phi = phi_cur;
                }
                skipCount = 0;
                normalizeVecs(slaveTrajBw, dimension, essDim, eps);
            }
            for (int i = 0; i < essDim; i++) {
                for (int j = 0; j < dimension; j++)
                    slaveTrajBw[i * dimension + j] += (traj_dots+stepCountPoints*dimension-dimension)[j];
            }
            skipCount++;

        }
        for (int i = 0; i < essDim; i++) {
            lyapExpBw[i] = lyapExpBw[i] / calcTime;
        }
    }
    else {
        for (int t = 0; t < calcTime; t += 1) {
            skipCount++;
            StepMAP(mainTrajectory, diffFunc, params, memI->arg, dimension);

            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajectories, dimension, expsNum, memI->projSum);
                ortVecs(essCur, dimension, essDim, memI->projSum);

                for (int i = 0; i < expsNum; i++)
                    lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / 1);

                normalizeVecs(slaveTrajectories, dimension, expsNum, 1);
                normalizeVecs(essCur, dimension, essDim, 1);
                memcpy(ncu_s + stepCountVectors * dimension * essDim, essCur, sizeof(double) * dimension * essDim);
                stepCountVectors++;
                skipCount = 0;
            }
            memcpy(traj_dots + stepCountPoints * dimension, mainTrajectory, sizeof(double) * dimension);
            stepCountPoints++;
        }

        for (int i = 0; i < expsNum; i++) {
            lyapExp[i] = lyapExp[i] / calcTime;
        }

        memcpy(main_afterForward, mainTrajectory, sizeof(double) * dimension);
        for (int t = 0; t < timeSkipClP; t += 1) {
            StepMAP(main_afterForward, diffFunc, params, memI->arg, dimension);
            memcpy(traj_dots + stepCountPoints * dimension, main_afterForward, sizeof(double) * dimension);
            stepCountPoints++;
        }
        // bw path
        double slaveTrajBw[MAX_DIMENSION * MAX_DIMENSION]; // exp_num * max_dimension
        memset(slaveTrajBw, 0., dimension * essDim * sizeof(double));
        for (int i = 0; i < essDim; ++i) {
            slaveTrajBw[i * dimension + i] += 1;
        }

        skipCount = 0;
        for (int t = 0; t < timeSkipClP; t += 1) {
            --stepCountPoints;
            if (skipCount == addSkipSteps) {
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                skipCount = 0;
            }
            skipCount++;
        }
        for (int t = 0; t < calcTime; t += 1) {
            stepCountPoints--;
            if (skipCount == addSkipSteps) {
                stepCountVectors--;
                ortVecs(slaveTrajBw, dimension, essDim, memI->projSum);
                normalizeVecs(slaveTrajBw, dimension, essDim, 1);
                for (int j = 0; j < essDim; j++)
                    for (int i = 0; i < essDim; i++)
                        memM->dotProductsMatrix[j*essDim + i] = dotProduct(ncu_s + stepCountVectors * dimension * essDim + j * dimension, slaveTrajBw + i * dimension, dimension);
                memcpy(memM->dotProductsMatrixT, memM->dotProductsMatrix, essDim * essDim * sizeof(double));
                matrixTranspose(memM->dotProductsMatrixT, essDim, essDim);
                matrixMult(memM->A, memM->dotProductsMatrixT, essDim, essDim, memM->dotProductsMatrix, essDim, essDim);
                double angle = pow(findEigenValues(memM->A, memM->R, memM->Q, memI->projSum, essDim, essDim, eigenIteration, memM->eigenValues), 0.5);
                double phi_cur;
                phi_cur = std::abs(3.14159265358979323846264338327950288 / 2. - acos(angle));
                if (phi_cur < *phi) {
                    *phi = phi_cur;
                }
                skipCount = 0;
            }
            skipCount++;
        }
    }
    return 0;
}
