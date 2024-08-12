#include <cmath>


void Shimizu_3D_flow(const double* state, double* res, const double* params) {
    res[0] = state[1] ;
    res[1] = - params[0] * state[1] - state[0] * state[2] + state[0] ;
    res[2] = - params[1] * state[2] + pow( state[0] , 2.0 ) ;
}

void Shimizu_3D_flow_var(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = state[1] *  1  ;
   res[1] = state[0] * ( 1 - stateOld[2] ) + state[1] * ( - params[0] ) + state[2] * ( - stateOld[0] ) ;
   res[2] = state[0] * ( 2.0 * pow( stateOld[0] , 1.0 ) ) + state[2] * ( - params[1] ) ;
}

void Shimizu_3D_flow_var_trans(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = state[1] *  1  ;
   res[1] = state[0] * ( 1 - stateOld[2] ) + state[1] * ( - params[0] ) + state[2] * ( - stateOld[0] ) ;
   res[2] = state[0] * ( 2.0 * pow( stateOld[0] , 1.0 ) ) + state[2] * ( - params[1] ) ;
}

void Shimizu_3D_flow_var_rev(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = -1 * ( state[0] * ( 0 ) +state[1] * ( 1 ) +state[2] * ( 0 )  );
   res[1] = -1 * ( state[0] * ( 1 - stateOld[2] ) +state[1] * ( - params[0] ) +state[2] * ( - stateOld[0] )  );
   res[2] = -1 * ( state[0] * ( 2.0 * pow( stateOld[0] , 1.0 ) ) +state[1] * ( 0 ) +state[2] * ( - params[1] )  );
}

void Shimizu_3D_flow_var_trans_rev(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = -1 * (state[1] * ( 1 - stateOld[2] ) + state[2] * ( 2.0 * pow( stateOld[0] , 1.0 ) ) );
   res[1] = -1 * (state[0] *  1  + state[1] * ( - params[0] ) );
   res[2] = -1 * (state[1] * ( - stateOld[0] ) + state[2] * ( - params[1] ) );
}
