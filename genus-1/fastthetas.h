#ifndef FAST_THETAS_H_
#define FAST_THETAS_H_

/* fastthetas.h -- headers for fastthetas.c
 *
 */


#ifdef __cplusplus
extern "C" {
#endif


// Cutoff point = 9000 digits
#ifndef LOW_PREC_THETA
#define LOW_PREC_THETA 30000
#endif

int distanceSmallerThan(mpc_t arg1, mpc_t arg2, mpfr_t bound);
void Pow2PowN(mpc_t r, mpc_t x, unsigned long int n);

void AGMPrime(mpc_t a[4]);
void Finfty4args(mpc_t lambda, mpc_t mu, mpc_t a[4]);
void Finfty2args(mpc_t lambda, mpc_t mu, mpc_t quo1, mpc_t quo2);
void WholeFunction(mpc_t lambda, mpc_t mu, mpc_t quo1, mpc_t quo2);
void NewtonStep(mpc_t quo2p[2], mpc_t quop[2], mpc_t lambda, mpc_t mu);
void MathfrakF(mpc_t thetas[4], mpc_t z, mpc_t tau);

void FastThetas(mpc_t thetas[6], mpc_t z, mpc_t tau);




#ifdef __cplusplus
}
#endif


#endif
