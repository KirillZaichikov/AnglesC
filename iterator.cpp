void StepMAP(double* val, void(*diffFunc)(const double*, double*, const double*), const double* params,
             double* arg, const int dimension){
    diffFunc(val, arg, params);
    for (int ii = 0; ii < dimension; ii++)
        val[ii] = arg[ii];
}

void StepMAPVAR(double* val, void(*diffFuncVar)(const double*, double*, const double*, const double*),
                const double *params, double* arg, const int dimension, const double* main){
    diffFuncVar(val, arg, params, main);
    for (int ii = 0; ii < dimension; ii++)
        val[ii] = arg[ii];
}