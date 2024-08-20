// #include "functionsLaser.cpp"
// #include "functions.cpp"
#include "functionsShimizu.cpp" // Нужно импортировать файл с функциями системы
#include "consts.h"
#include "Node.cpp"
#include "utils.h"
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

// Дискретизация плоскости параметров
#define DISCR 1 

int main(){
    printf("***PROGRAM START***\n");
    // TODO: адекватный ввод функций
    void (*diffFunc)(const double*, double*, const double*) {Shimizu_3D_flow};
    void (*diffFuncRev)(const double*, double*, const double*) {Shimizu_3D_flow_rev};
    void (*diffFuncVar)(const double*, double*, const double*, const double*) {Shimizu_3D_flow_var};
    void (*diffFuncVarTrans)(const double*, double*, const double*, const double*) {Shimizu_3D_flow_var_trans};
    void (*diffFuncVarRev)(const double*, double*, const double*, const double*) {Shimizu_3D_flow_var_rev};
    void (*diffFuncVarTransRev)(const double*, double*, const double*, const double*) {Shimizu_3D_flow_var_trans_rev};
    // void (*diffFunc)(const double*, double*, const double*) {laser_6D_flow};
    // void (*diffFuncVar)(const double*, double*, const double*, const double*) {laser_6D_flow_var};
    // void (*diffFuncVarTrans)(const double*, double*, const double*, const double*) {laser_6D_flow_var_trans};
    // void (*diffFuncVarRev)(const double*, double*, const double*, const double*) {laser_6D_flow_var_rev};
    // void (*diffFuncVarTransRev)(const double*, double*, const double*, const double*) {laser_6D_flow_var_trans_rev};
    
    // Параметры можно вводить как точкой так и плоскостью
    int paramsDim = 2;
    double params[] = {0.75, 0.5};
    //double params[] = {1, -1, 1.8, 0.8, 0.5};
    // double params[] = {0, 0, 50.0, -1.0, 0.0, 1.15};
    double param1[DISCR];
    double param2[DISCR];
    // for (int i=0;i<DISCR;i++){
    //     param1[i] = 4.666 + (4.710-4.666)*i/DISCR;
    //     param2[i] = 0.117 + (0.152-0.117)*i/DISCR;
    //     std::cout << param1[i];
    // }
    param1[0] = 0.75;
    param2[0] = 0.5;

    ofstream fout;
    fout.open("File.txt");
    
    bool bySlaveTraj = 1; 
    bool type = 1; // 1 - flow

    int samplesCounter = 0;
    
    // TODO: Нужно параметры для каждой системы заполнять автоматически
    int dimension = 3;
    int expsNum = 3;
    int essDim = 1;
    int samplesNum = 1;
    double eps = 0.00001; // В вариациях не имеет значения
    double step = 0.01;
    double timeSkip = 100;
    double timeSkipClP = 50;
    double calcTime = 2000;
    double addSkipSteps = 1;
    double inits[] = {0.0000001, 0, 0, 0, 0};

    double main[MAX_DIMENSION];
    double lExp[MAX_DIMENSION];
    double lExp1[MAX_DIMENSION];
    double lExpBw[MAX_DIMENSION];
    double lExpBw1[MAX_DIMENSION];
    double lExpT[MAX_DIMENSION];
    double lExpT1[MAX_DIMENSION];
    double sub[MAX_DIMENSION * MAX_DIMENSION];
    double ess[MAX_DIMENSION * MAX_DIMENSION];
    double minPhi = 300;
    MemoryForIntegrator integratorParams{};
    MemoryForMatrix masForMatrix{};

    int eigenIteration;
    if ((essDim == 1) || (dimension-essDim == 1)){
        eigenIteration = 1;
    } else {
        eigenIteration = 50;
    }

    // Динмаические массивы для хранения точек и векторов
    size_t traj_dots_cnt = (size_t) ((calcTime / samplesNum + timeSkipClP) / step);
    double* traj_s = (double*)malloc(traj_dots_cnt * dimension * sizeof(double));
    size_t ncu_dots_cnt = (size_t) (calcTime / samplesNum / step);
    double* ncu_s;

    // Определяем длину массива векторов и уведомляем пользователя о том какой алгоритм выбран
    if (bySlaveTraj){
        if ((dimension / 2 <= essDim) && (essDim!=1)){
            std::cout << "Slave Trajectory calc by CU vectors" << std::endl;
            ncu_s = (double*)malloc(ncu_dots_cnt * dimension * (dimension - essDim) * sizeof(double));
        } else {
            std::cout << "Slave Trajectory calc by SS vectors" << std::endl;
            ncu_s = (double*)malloc(ncu_dots_cnt * dimension * essDim * sizeof(double));
        }
    }
    else{ 
        if ((dimension / 2 <= essDim) && (essDim !=1)){
            std::cout << "VAR calc by CU vectors" << std::endl;
            ncu_s = (double*)malloc(ncu_dots_cnt * dimension * (dimension - essDim) * sizeof(double));
        } else {
            std::cout << "VAR calc by SS vectors" << std::endl;
            ncu_s = (double*)malloc(ncu_dots_cnt * dimension * essDim * sizeof(double));
        }
    }

    // Здесь должно быть два цикла for, чтобы пробегать по параметрам, пока убрал, чтобы тестить было удобнее
    printf("%lf %lf\n", param1[0], param2[0]);
    memset(lExp1, 0., MAX_DIMENSION * sizeof(double));
    memset(lExpBw1, 0., MAX_DIMENSION * sizeof(double));
    minPhi=300;
    for (int d = 0; d < dimension; ++d) {
        main[d] = inits[d];
    }

    double paramsAll[] = {param1[0], param2[0]};
    if (bySlaveTraj) {
        if ((dimension / 2 <= essDim) && (essDim !=1)) {
            printf("***START WITH SLAVE SS***\n");
            cudaAngleSlaveNodeWarm(main, sub, ess, dimension, diffFunc, paramsAll, eps, expsNum,
                            step, timeSkip, timeSkipClP, type, essDim, &integratorParams, addSkipSteps);
            printf("End Warm\n");
            for (int s = 0; s < samplesNum; s++) {
                auto res = cudaAngleSlaveNodeCalculate(main, sub, ess, dimension,
                    diffFunc, diffFuncRev, paramsAll, eps, lExp, lExpBw, &minPhi,
                    expsNum, step, timeSkip, timeSkipClP, calcTime / samplesNum, addSkipSteps,
                    traj_s, ncu_s, type, essDim, &integratorParams, &masForMatrix, eigenIteration);
                for (int i =0; i<expsNum; i++){
                    lExp1[i] += lExp[i];
                    lExpBw1[i] += lExpBw[i];
                }
                samplesCounter++;
                if (res == -1)
                    break;
            }
        } else {
            
        }
    } else {
        if ((dimension / 2 <= essDim) && (essDim !=1)) {
            printf("***START WITH CU BASIS***\n");
            cudaAngleCUNodeWarm(main, sub, ess, dimension, diffFunc, diffFuncVar, paramsAll, eps, expsNum,
                            step, timeSkip, timeSkipClP, addSkipSteps, type, essDim, &integratorParams);
            for (int s = 0; s < samplesNum; s++) {
                auto res = cudaAngleCUNodeCalculate(main, sub, ess, dimension,
                    diffFunc, diffFuncVar, diffFuncVarTrans,
                    paramsAll, eps, lExp, lExpBw, &minPhi,
                    expsNum, step, timeSkip, timeSkipClP, calcTime / samplesNum, addSkipSteps,
                    traj_s, ncu_s, type, essDim, &integratorParams, &masForMatrix, eigenIteration);
                for (int i =0; i<expsNum; i++){
                    lExp1[i] += lExp[i];
                    lExpBw1[i] += lExpBw[i];
                }
                samplesCounter++;
                if (res == -1)
                    break;
            }
        } else {
            printf("***START WITH SS BASIS***\n");
            cudaAngleSSNodeWarm(main, sub, ess, dimension, diffFunc, diffFuncVar, diffFuncVarTransRev, paramsAll, eps, expsNum,
                            step, timeSkip, timeSkipClP, type, essDim, &integratorParams, addSkipSteps);
            for (int s = 0; s < samplesNum; s++) {
                auto res = cudaAngleSSNodeCalculate(main, sub, ess, dimension,
                    diffFunc, diffFuncVar, diffFuncVarRev, diffFuncVarTransRev,
                    paramsAll, eps, lExp, lExpBw, &minPhi,
                    expsNum, step, timeSkip, timeSkipClP, calcTime / samplesNum, addSkipSteps,
                    traj_s, ncu_s, type, essDim, &integratorParams, &masForMatrix, eigenIteration);
                for (int i =0; i<expsNum; i++){
                    lExp1[i] += lExp[i];
                    lExpBw1[i] += lExpBw[i];
                    lExpT1[i] +=lExpT[i];
                }
                samplesCounter++;
                if (res == -1)
                    break;
            }
        }
    }

    fout<<param1[0]<<" "<<param2[0]<<" "<<lExp1[0]/samplesCounter<<" "<< minPhi <<", ";

    for (int i =0; i<expsNum; i++)
        std::cout << lExp1[i]/(samplesCounter) << " ";
    std::cout  << std::endl;
    for (int i =0; i<essDim; i++)
        std::cout << lExpBw1[i]/(samplesCounter) << " ";
    std::cout << std::endl;
    for (int i =0; i<expsNum; i++)
        std::cout << lExpT1[i]/(samplesCounter) << " ";
    std::cout << std::endl;
    std::cout << "Angle:" << minPhi;
    return 0;
}

void memSetZero(struct MemoryForIntegrator intMem, struct MemoryForMatrix matMem) {
    /*После работы алгоритма вдоль одного куска зануляет переменные*/
    memset(matMem.dotProductsMatrix, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(matMem.dotProductsMatrixT, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(matMem.A, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(matMem.R, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(matMem.Q, 0., MAX_DIMENSION * MAX_DIMENSION * sizeof(double));
    memset(matMem.eigenValues, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.k1, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.k2, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.k3, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.k4, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.k5, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.k6, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.k7, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.k8, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.arg, 0., MAX_DIMENSION * sizeof(double));
    memset(intMem.projSum, 0., MAX_DIMENSION * sizeof(double));
}