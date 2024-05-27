#include <cmath>


void LorenzDoub_6D_flow(const double* state, double* res, const double* params) {
    res[0] = - params[4] * ( state[0] - state[3] ) - params[0] * ( state[0] - state[1] ) ;
    res[1] = state[0] * ( params[1] - state[2] ) - state[1] ;
    res[2] = - params[2] * state[2] + state[0] * state[1] ;
    res[3] = params[4] * ( state[0] - state[3] ) - params[0] * ( state[3] - state[4] ) ;
    res[4] = - state[3] * state[5] + state[3] * ( params[3] + params[1] ) - state[4] ;
    res[5] = - params[2] * state[5] + state[3] * state[4] ;
}

void LorenzDoub_6D_flow_var(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = state[0] * ( - params[4] - params[0] ) + state[1] * ( params[0] ) + state[3] * ( params[4] ) ;
   res[1] = state[0] * ( params[1] - stateOld[2] ) + state[1] * ( - 1 ) + state[2] * ( - stateOld[0] ) ;
   res[2] = state[0] * ( stateOld[1] ) + state[1] * ( stateOld[0] ) + state[2] * ( - params[2] ) ;
   res[3] = state[0] * ( params[4] ) + state[3] * ( - params[4] - params[0] ) + state[4] * ( params[0] ) ;
   res[4] = state[3] * ( params[3] + params[1] - stateOld[5] ) + state[4] * ( - 1 ) + state[5] * ( - stateOld[3] ) ;
   res[5] = state[3] * ( stateOld[4] ) + state[4] * ( stateOld[3] ) + state[5] * ( - params[2] ) ;
}

void LorenzDoub_6D_flow_var_trans(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = state[0] * ( - params[4] - params[0] ) + state[1] * ( params[1] - stateOld[2] ) + state[2] * ( stateOld[1] ) + state[3] * ( params[4] ) ;
   res[1] = state[0] * ( params[0] ) + state[1] * ( - 1 ) + state[2] * ( stateOld[0] ) ;
   res[2] = state[1] * ( - stateOld[0] ) + state[2] * ( - params[2] ) ;
   res[3] = state[0] * ( params[4] ) + state[3] * ( - params[4] - params[0] ) + state[4] * ( params[3] + params[1] - stateOld[5] ) + state[5] * ( stateOld[4] ) ;
   res[4] = state[3] * ( params[0] ) + state[4] * ( - 1 ) + state[5] * ( stateOld[3] ) ;
   res[5] = state[4] * ( - stateOld[3] ) + state[5] * ( - params[2] ) ;
}

void LorenzDoub_6D_flow_var_rev(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = -1 * ( state[0] * ( - params[4] - params[0] ) +state[1] * ( params[0] ) +state[2] * ( 0 ) +state[3] * ( params[4] ) +state[4] * ( 0 ) +state[5] * ( 0 )  );
   res[1] = -1 * ( state[0] * ( params[1] - stateOld[2] ) +state[1] * ( - 1 ) +state[2] * ( - stateOld[0] ) +state[3] * ( 0 ) +state[4] * ( 0 ) +state[5] * ( 0 )  );
   res[2] = -1 * ( state[0] * ( stateOld[1] ) +state[1] * ( stateOld[0] ) +state[2] * ( - params[2] ) +state[3] * ( 0 ) +state[4] * ( 0 ) +state[5] * ( 0 )  );
   res[3] = -1 * ( state[0] * ( params[4] ) +state[1] * ( 0 ) +state[2] * ( 0 ) +state[3] * ( - params[4] - params[0] ) +state[4] * ( params[0] ) +state[5] * ( 0 )  );
   res[4] = -1 * ( state[0] * ( 0 ) +state[1] * ( 0 ) +state[2] * ( 0 ) +state[3] * ( params[3] + params[1] - stateOld[5] ) +state[4] * ( - 1 ) +state[5] * ( - stateOld[3] )  );
   res[5] = -1 * ( state[0] * ( 0 ) +state[1] * ( 0 ) +state[2] * ( 0 ) +state[3] * ( stateOld[4] ) +state[4] * ( stateOld[3] ) +state[5] * ( - params[2] )  );
}

void LorenzDoub_6D_flow_var_trans_rev(const double* state, double* res, const double* params, const double* stateOld) {
   res[0] = -1 * (state[0] * ( - params[4] - params[0] ) + state[1] * ( params[1] - stateOld[2] ) + state[2] * ( stateOld[1] ) + state[3] * ( params[4] ) );
   res[1] = -1 * (state[0] * ( params[0] ) + state[1] * ( - 1 ) + state[2] * ( stateOld[0] ) );
   res[2] = -1 * (state[1] * ( - stateOld[0] ) + state[2] * ( - params[2] ) );
   res[3] = -1 * (state[0] * ( params[4] ) + state[3] * ( - params[4] - params[0] ) + state[4] * ( params[3] + params[1] - stateOld[5] ) + state[5] * ( stateOld[4] ) );
   res[4] = -1 * (state[3] * ( params[0] ) + state[4] * ( - 1 ) + state[5] * ( stateOld[3] ) );
   res[5] = -1 * (state[4] * ( - stateOld[3] ) + state[5] * ( - params[2] ) );
}