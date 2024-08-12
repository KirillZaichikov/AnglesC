#include <cmath>


void Gonchenko_3D_map(const double* state, double* res, const double* params) {
    res[0] = state[1] ;
    res[1] = ( - params[4] * state[0] + params[3] * state[1] * ( params[0] * state[2] + params[2] ) ) / params[3] ;
    res[2] = params[1] * pow ( state[1] , 2.0 ) + params[3] * state[2] ;
}

void Gonchenko_3D_map_var(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = state[1] *  1  ;
   res[1] = state[0] * ( - params[4] / params[3] ) + state[1] * ( params[0] * stateOld[2] + params[2] ) + state[2] * ( params[0] * stateOld[1] ) ;
   res[2] = state[1] * ( 2.0 * params[1] * pow ( stateOld[1] , 1.0 ) ) + state[2] * ( params[3] ) ;
}

void Gonchenko_3D_map_var_trans(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = state[1] * ( - params[4] / params[3] ) ;
   res[1] = state[0] *  1  + state[1] * ( params[0] * stateOld[2] + params[2] ) + state[2] * ( 2.0 * params[1] * pow ( stateOld[1] , 1.0 ) ) ;
   res[2] = state[1] * ( params[0] * stateOld[1] ) + state[2] * ( params[3] ) ;
}

void Gonchenko_3D_map_var_rev(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = state[0] * ( ( - 2.0 * params[0] * params[1] * pow ( stateOld[0] , 2.0 ) - 1.0 * params[0] * ( params[1] * pow ( stateOld[0] , 2.0 ) - stateOld[2] ) + params[2] * params[3] ) / params[4] ) + state[1] * ( - params[3] / params[4] ) + state[2] * ( 1.0 * params[0] * stateOld[0] / params[4] ) ;
   res[1] = state[0] *  1  ;
   res[2] = state[0] * ( - 2.0 * params[1] * pow ( stateOld[0] , 1.0 ) / params[3] ) + state[2] * ( 1.0 / params[3] ) ;
}

void Gonchenko_3D_map_var_trans_rev(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = state[0] * ( ( - 2.0 * params[0] * params[1] * pow ( stateOld[0] , 2.0 ) - 1.0 * params[0] * ( params[1] * pow ( stateOld[0] , 2.0 ) - stateOld[2] ) + params[2] * params[3] ) / params[4] ) + state[1] *  1  + state[2] * ( - 2.0 * params[1] * pow ( stateOld[0] , 1.0 ) / params[3] ) ;
   res[1] = state[0] * ( - params[3] / params[4] ) ;
   res[2] = state[0] * ( 1.0 * params[0] * stateOld[0] / params[4] ) + state[2] * ( 1.0 / params[3] ) ;
}