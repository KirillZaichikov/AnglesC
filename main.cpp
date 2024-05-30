// #include "functionsLaser.cpp"
#include "functions.cpp"
#include "consts.h"
#include "Node.cpp"
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#define DISCR 200
using namespace std;

int main(){
    printf("***PROGRAM START***\n");
    void (*diffFunc)(const double*, double*, const double*) {LorenzDoub_6D_flow};
    void (*diffFuncVar)(const double*, double*, const double*, const double*) {LorenzDoub_6D_flow_var};
    void (*diffFuncVarTrans)(const double*, double*, const double*, const double*) {LorenzDoub_6D_flow_var_trans};
    void (*diffFuncVarRev)(const double*, double*, const double*, const double*) {LorenzDoub_6D_flow_var_rev};
    void (*diffFuncVarTransRev)(const double*, double*, const double*, const double*) {LorenzDoub_6D_flow_var_trans_rev};
    // void (*diffFunc)(const double*, double*, const double*) {laser_6D_flow};
    // void (*diffFuncVar)(const double*, double*, const double*, const double*) {laser_6D_flow_var};
    // void (*diffFuncVarTrans)(const double*, double*, const double*, const double*) {laser_6D_flow_var_trans};
    // void (*diffFuncVarRev)(const double*, double*, const double*, const double*) {laser_6D_flow_var_rev};
    // void (*diffFuncVarTransRev)(const double*, double*, const double*, const double*) {laser_6D_flow_var_trans_rev};
    
    int paramsDim = 6;
    double params[] = {9, 25, 2.666666666666, 0, 0.1};
    // double params[] = {0, 0, 50.0, -1.0, 0.0, 1.15};
    double param1[DISCR];
    double param2[DISCR];
    // for (int i=0;i<DISCR;i++){
    //     param1[i] = 4.666 + (4.710-4.666)*i/DISCR;
    //     param2[i] = 0.117 + (0.152-0.117)*i/DISCR;
    //     std::cout << param1[i];
    // }
    param1[0] = 9;
    param2[0] = 25;
    ofstream fout;
    fout.open("File.txt");
    
    int dimension = 6;
    int samplesNum = 1;
    int samplesCounter = 0;
    int expsNum = 6;
    int essDim = 2;
    double eps = 1;
    double step = 0.001;
    double timeSkip = 100;
    double timeSkipClP = 50;
    double calcTime = 1000;
    double addSkipSteps = 1;
    double inits[] = {0.0000001, 0, 0, 0, 0, 0};
    double main[MAX_DIMENSION];
    double lExp[MAX_DIMENSION];
    double lExpBw[MAX_DIMENSION];
    double lExp1[MAX_DIMENSION];
    double lExpBw1[MAX_DIMENSION];
    double lExpT[MAX_DIMENSION];
    double lExpT1[MAX_DIMENSION];
    double sub[MAX_DIMENSION * MAX_DIMENSION];
    double ess[MAX_DIMENSION * MAX_DIMENSION];
    double minPhi = 300;

    size_t traj_dots_cnt = (size_t) ((calcTime / samplesNum + timeSkipClP) / step);
    double* traj_s = (double*)malloc(traj_dots_cnt * dimension * sizeof(double));
    size_t ncu_dots_cnt = (size_t) (calcTime / samplesNum / step);
    double* ncu_s;
    if (dimension / 2 <= essDim){
        std::cout << "CALC BY CU VECTORS" << std::endl;
        ncu_s = (double*)malloc(ncu_dots_cnt * dimension * (dimension - essDim) * sizeof(double));
    } else {
        std::cout << "Calc By SS vectors" << std::endl;
        ncu_s = (double*)malloc(ncu_dots_cnt * dimension * essDim * sizeof(double));
    }

    printf("%lf %lf\n", param1[0], param2[0]);
    memset(lExp1, 0., MAX_DIMENSION * sizeof(double));
    memset(lExpBw1, 0., MAX_DIMENSION * sizeof(double));
    minPhi=300;
    for (int d = 0; d < dimension; ++d) {
        main[d] = inits[d];
    }
    double paramsAll[] = {param1[0], param2[0], params[2], params[3], params[4], params[5]};
    if (dimension / 2 <= essDim) {
        cudaAngleNode2Warm(main, sub, ess, dimension, diffFunc, diffFuncVar, paramsAll, eps, expsNum,
                        step, timeSkip, timeSkipClP, essDim);
        for (int s = 0; s < samplesNum; s++) {
            auto res = cudaAngleNode2Calculate(main, sub, ess, dimension,
                diffFunc, diffFuncVar, diffFuncVarTrans,
                paramsAll, eps, lExp, lExpBw, &minPhi,
                expsNum, step, timeSkip, timeSkipClP, calcTime / samplesNum, addSkipSteps,
                traj_s, ncu_s, essDim);
            for (int i =0; i<expsNum; i++){
                lExp1[i] += lExp[i];
                lExpBw1[i] += lExpBw[i];
            }
            samplesCounter++;
            if (res == -1)
                break;
        }
    } else {
        printf("***START WARM***\n");
        cudaAngleNodeWarm(main, sub, ess, dimension, diffFunc, diffFuncVar, diffFuncVarTransRev, params, eps, expsNum,
                        step, timeSkip, timeSkipClP, essDim);
        printf("***END WARM***\n");
        for (int s = 0; s < samplesNum; s++) {
            auto res = cudaAngleNodeCalculate(main, sub, ess, dimension,
                diffFunc, diffFuncVar, diffFuncVarRev, diffFuncVarTransRev,
                params, eps, lExp, lExpBw, lExpT, &minPhi,
                expsNum, step, timeSkip, timeSkipClP, calcTime / samplesNum, addSkipSteps,
                traj_s, ncu_s, essDim);
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