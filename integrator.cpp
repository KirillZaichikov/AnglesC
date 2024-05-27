void dverkStep(double* val, const int dimension, void(*diffFunc)(const double*, double*, const double*), const double* params, double step,
    double* arg, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double* k8) {

    diffFunc(val, k1, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step / 6. * k1[j];
    }
    diffFunc(arg, k2, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (4. / 75. * k1[j] + 16. / 75. * k2[j]);
    }
    diffFunc(arg, k3, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (5. / 6. * k1[j] - 8. / 3. * k2[j] + 5. / 2. * k3[j]);
    }
    diffFunc(arg, k4, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (-165. / 64. * k1[j] + 55. / 6. * k2[j] - 425. / 64. * k3[j] + 85. / 96. * k4[j]);
    }
    diffFunc(arg, k5, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (12. / 5. * k1[j] - 8. * k2[j] + 4015. / 612. * k3[j] - 11. / 36. * k4[j] + 88. / 255. * k5[j]);
    }
    diffFunc(arg, k6, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (-8263. / 15000. * k1[j] + 124. / 75. * k2[j] - 643. / 680. * k3[j] - 81. / 250. * k4[j] + 2484. / 10625. * k5[j]);
    }
    diffFunc(arg, k7, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (3501. / 1720. * k1[j] - 300. / 43. * k2[j] + 297275. / 52632. * k3[j] - 319. / 2322. * k4[j] + 24068. / 84065. * k5[j] + 3850. / 26703. * k7[j]);
    }
    diffFunc(arg, k8, params);
    for (int j = 0; j < dimension; j++) {
        val[j] = val[j] + step * (3. / 40. * k1[j] + 875. / 2244. * k3[j] + 23. / 72. * k4[j] + 264. / 1955. * k5[j] + 125. / 11592. * k7[j] + 43. / 616. * k8[j]);
    }
}


void dverkStepVarMat(double* val, const int dimension, void(*diffFuncVar)(const double*, double*, const double*, const double*), const double* params, double step,
    double* arg, const double* mainTraj, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double* k8) {

    diffFuncVar(val, k1, params, mainTraj);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step / 6. * k1[j];
    }
    diffFuncVar(arg, k2, params, mainTraj);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (4. / 75. * k1[j] + 16. / 75. * k2[j]);
    }
    diffFuncVar(arg, k3, params, mainTraj);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (5. / 6. * k1[j] - 8. / 3. * k2[j] + 5. / 2. * k3[j]);
    }
    diffFuncVar(arg, k4, params, mainTraj);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (-165. / 64. * k1[j] + 55. / 6. * k2[j] - 425. / 64. * k3[j] + 85. / 96. * k4[j]);
    }
    diffFuncVar(arg, k5, params, mainTraj);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (12. / 5. * k1[j] - 8. * k2[j] + 4015. / 612. * k3[j] - 11. / 36. * k4[j] + 88. / 255. * k5[j]);
    }
    diffFuncVar(arg, k6, params, mainTraj);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (-8263. / 15000. * k1[j] + 124. / 75. * k2[j] - 643. / 680. * k3[j] - 81. / 250. * k4[j] + 2484. / 10625. * k5[j]);
    }
    diffFuncVar(arg, k7, params, mainTraj);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (3501. / 1720. * k1[j] - 300. / 43. * k2[j] + 297275. / 52632. * k3[j] - 319. / 2322. * k4[j] + 24068. / 84065. * k5[j] + 3850. / 26703. * k7[j]);
    }
    diffFuncVar(arg, k8, params, mainTraj);
    for (int j = 0; j < dimension; j++) {
        val[j] = val[j] + step * (3. / 40. * k1[j] + 875. / 2244. * k3[j] + 23. / 72. * k4[j] + 264. / 1955. * k5[j] + 125. / 11592. * k7[j] + 43. / 616. * k8[j]);
    }
}