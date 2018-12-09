#ifndef NAIVE_THETAS_H_
#define NAIVE_THETAS_H_

/* naivethetas.h -- headers for naivethetas.c
 *
 */


#ifdef __cplusplus
extern "C" {
#endif




void NaiveTheta(mpc_t res, mpc_t z, mpc_t tau, unsigned long int p);
void InternalNaiveThetaFunctionAndThetaConstants_00_01(mpc_t res[4], mpc_t z, mpc_t tau, unsigned long int p);
void NaiveThetaFunctionAndThetaConstants_00_01(mpc_t res[4], mpc_t z, mpc_t tau, unsigned long int p);
void NaiveThetaFunctionAndThetaConstants_00_01_10(mpc_t res[6], mpc_t z, mpc_t tau, unsigned long int p);

void NaiveThetaConstants_00_01(mpc_t res[2], mpc_t tau, unsigned long int p);


#ifdef __cplusplus
}
#endif


#endif
