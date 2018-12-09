#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mpc.h>
#include <stdlib.h>
#include <sys/types.h>          /* for cputime */
#include <sys/resource.h>       /* for cputime */

#include "naivethetas.h"
#include "fastthetas.h"
#include "fastthetaconstants.h"

uint64_t microseconds()
{
    struct rusage res[1];
    getrusage(RUSAGE_SELF, res);
    uint64_t r;
    r = (uint64_t) res->ru_utime.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) res->ru_utime.tv_usec;
    return r;
}

double cputime()
{
    double milli = (microseconds() / (uint64_t) 1000);
    return milli / 1.0e3;
}

void usage(const char * progname) {
    fprintf(stderr, "Usage: %s [function] [precision]\n", progname);
    exit(1);
}


int main(int argc, char * argv[])
{
    int add_time = 0;
    double tt = 0;
    int i;

    if (argc < 2) usage(argv[0]);

    if (strcmp(argv[1], "time") == 0) {
        add_time = 1;
        argv++,argc--;
    }

    if (argc != 3) usage(argv[0]);

    int prec = atoi(argv[2]);

    if (prec == 0) usage(argv[0]);


    mpc_t agmprime[4]; for (i=0; i<4; i++) mpc_init2(agmprime[i], prec);
    mpc_t r; mpc_init2(r, prec);
    mpc_t x; mpc_init2(x, prec);
    unsigned long n;
    mpc_t lambda; mpc_init2(lambda, prec);
    mpc_t mu; mpc_init2(mu, prec);
    mpc_t argF[4]; for (i=0; i<4; i++) mpc_init2(argF[i], prec);
    mpc_t quo1; mpc_init2(quo1, prec);
    mpc_t quo2; mpc_init2(quo2, prec);
    mpc_t quop[2]; for (i=0; i<2; i++) mpc_init2(quop[i], prec);
    mpc_t quo2p[2]; for (i=0; i<2; i++) mpc_init2(quo2p[i], prec);
    mpc_t theta; mpc_init2(theta, prec);
    mpc_t z; mpc_init2(z, prec);
    mpc_t tau; mpc_init2(tau, prec);
    mpc_t theta00and01[4]; for (i=0; i<4; i++) mpc_init2(theta00and01[i], prec);
    mpc_t thetas[6]; for (i=0; i<6; i++) mpc_init2(thetas[i], prec);
    mpc_t thetaconstants[3]; for (i=0; i<3; i++) mpc_init2(thetaconstants[i], prec);
    mpc_t constants2[2]; for (i=0; i<2; i++) mpc_init2(constants2[i], prec);


    
    // These are the functions you're looking for
    if(strcmp(argv[1], "FastThetas") == 0) {
        size_t read = 0;
        mpc_inp_str(z, stdin, &read, 0, MPC_RNDNN);
        mpc_inp_str(tau, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        FastThetas(thetas, z, tau);
        tt += cputime();
        for(i = 0 ; i < 6 ; i++) {
            mpc_out_str(stdout, 10, 0, thetas[i], MPC_RNDNN);
            printf("\n");
        }
    } else if(strcmp(argv[1], "DupontThetaConstants") == 0) {
        size_t read = 0;
        mpc_inp_str(tau, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        DupontThetaConstants(thetaconstants, tau);
        tt += cputime();
        for(i = 0 ; i < 3 ; i++) {
            mpc_out_str(stdout, 10, 0, thetaconstants[i], MPC_RNDNN);
            printf("\n");
        }
    } else if(strcmp(argv[1], "NaiveTheta") == 0) {
        size_t read = 0;
        mpc_inp_str(z, stdin, &read, 0, MPC_RNDNN);
        mpc_inp_str(tau, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        NaiveTheta(theta, z, tau, mpfr_get_prec(mpc_realref(tau)));
        tt += cputime();
        mpc_out_str(stdout, 10, 0, theta, MPC_RNDNN);
        printf("\n");
    } else if(strcmp(argv[1], "NaiveThetaFunctionAndThetaConstants_00_01") == 0) {
        size_t read = 0;
        mpc_inp_str(z, stdin, &read, 0, MPC_RNDNN);
        mpc_inp_str(tau, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        NaiveThetaFunctionAndThetaConstants_00_01(theta00and01, z, tau, mpfr_get_prec(mpc_realref(tau))-10);
        tt += cputime();
        for(i = 0 ; i < 4 ; i++) {
            mpc_out_str(stdout, 10, 0, theta00and01[i], MPC_RNDNN);
            printf("\n");
        }
    } else if(strcmp(argv[1], "NaiveThetaFunctionAndThetaConstants_00_01_10") == 0) {
        size_t read = 0;
        mpc_inp_str(z, stdin, &read, 0, MPC_RNDNN);
        mpc_inp_str(tau, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        NaiveThetaFunctionAndThetaConstants_00_01_10(thetas, z, tau, mpfr_get_prec(mpc_realref(tau))-1);
        tt += cputime();
        for(i = 0 ; i < 6 ; i++) {
            mpc_out_str(stdout, 10, 0, thetas[i], MPC_RNDNN);
            printf("\n");
        }
    } else if(strcmp(argv[1], "NaiveThetaConstants_00_01") == 0) {
        size_t read = 0;
        mpc_inp_str(tau, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        NaiveThetaConstants_00_01(constants2, tau, mpfr_get_prec(mpc_realref(tau)));
        tt += cputime();
        for(i = 0 ; i < 2 ; i++) {
            mpc_out_str(stdout, 10, 0, constants2[i], MPC_RNDNN);
            printf("\n");
        }
        
        
    // Other functions
    } else if (strcmp(argv[1], "AGMPrime") == 0) {
        /* 4 arguments */
        size_t read = 0;
        for(i = 0 ; i < 4 ; i++) {
            mpc_inp_str(agmprime[i], stdin, &read, 0, MPC_RNDNN);
        }
        tt -= cputime();
        AGMPrime(agmprime);
        tt += cputime();
        for(i = 0 ; i < 4 ; i++) {
            mpc_out_str(stdout, 10, 0, agmprime[i], MPC_RNDNN);
            printf("\n");
        }
    } else if (strcmp(argv[1], "Pow2PowN") == 0) {
        /* x and n */
        size_t read = 0;
        mpc_inp_str(x, stdin, &read, 0, MPC_RNDNN);
        scanf("%lu", &n);
        tt -= cputime();
        Pow2PowN(r, x, n);
        tt += cputime();
        mpc_out_str(stdout, 10, 0, r, MPC_RNDNN);
        printf("\n");
    } else if(strcmp(argv[1], "Finfty4args") == 0) {
        /* array of 4 arguments */
        size_t read = 0;
        for(i = 0 ; i < 4 ; i++) {
            mpc_inp_str(argF[i], stdin, &read, 0, MPC_RNDNN);
        }
        tt -= cputime();
        Finfty4args(lambda, mu, argF);
        tt += cputime();
        mpc_out_str(stdout, 10, 0, lambda, MPC_RNDNN);		printf("\n");
        mpc_out_str(stdout, 10, 0, mu, MPC_RNDNN);
        printf("\n");
    } else if(strcmp(argv[1], "Finfty2args") == 0) {
        /* 2 quotients */
        size_t read = 0;
        mpc_inp_str(quo1, stdin, &read, 0, MPC_RNDNN);
        mpc_inp_str(quo2, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        Finfty2args(lambda, mu, quo1, quo2);
        tt += cputime();
        mpc_out_str(stdout, 10, 0, lambda, MPC_RNDNN);		printf("\n");
        mpc_out_str(stdout, 10, 0, mu, MPC_RNDNN);
        printf("\n");
    } else if(strcmp(argv[1], "WholeFunction") == 0) {
        size_t read = 0;
        mpc_inp_str(quo1, stdin, &read, 0, MPC_RNDNN);
        mpc_inp_str(quo2, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        WholeFunction(lambda, mu, quo1, quo2);
        tt += cputime();
        mpc_out_str(stdout, 10, 0, lambda, MPC_RNDNN);		printf("\n");
        mpc_out_str(stdout, 10, 0, mu, MPC_RNDNN);
        printf("\n");
    } else if(strcmp(argv[1], "NewtonStep") == 0) {
        size_t read = 0;
        for(i = 0 ; i < 2 ; i++) {
            mpc_inp_str(quop[i], stdin, &read, 0, MPC_RNDNN);
        }
        mpc_inp_str(lambda, stdin, &read, 0, MPC_RNDNN);
        mpc_inp_str(mu, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        NewtonStep(quo2p, quop, lambda, mu);
        tt += cputime();
        for(i = 0 ; i < 2 ; i++) {
            mpc_out_str(stdout, 10, 0, quo2p[i], MPC_RNDNN);
            printf("\n");
        }
    } else if(strcmp(argv[1], "InternalNaiveThetaFunctionAndThetaConstants_00_01") == 0) {
        size_t read = 0;
        mpc_inp_str(z, stdin, &read, 0, MPC_RNDNN);
        mpc_inp_str(tau, stdin, &read, 0, MPC_RNDNN);
        tt -= cputime();
        InternalNaiveThetaFunctionAndThetaConstants_00_01(theta00and01, z, tau, mpfr_get_prec(mpc_realref(tau)));
        tt += cputime();
        for(i = 0 ; i < 4 ; i++) {
            mpc_out_str(stdout, 10, 0, theta00and01[i], MPC_RNDNN);
            printf("\n");
        }
    } else {
        usage(argv[0]);
    }

    if (add_time) {
        //printf("%1.7f\n", tt);
      mpc_set_d(r, tt, MPC_RNDNN);
      mpc_out_str(stdout, 10, 0, r, MPC_RNDNN);
      printf("\n");
    }

    for (i=0; i<4; i++) mpc_clear(agmprime[i]);
    mpc_clear(r);
    mpc_clear(x);
    mpc_clear(lambda);
    mpc_clear(mu);
    for (i=0; i<4; i++) mpc_clear(argF[i]);
    mpc_clear(quo1);
    mpc_clear(quo2);
    for (i=0; i<2; i++) mpc_clear(quop[i]);
    for (i=0; i<2; i++) mpc_clear(quo2p[i]);
    mpc_clear(theta);
    mpc_clear(z);
    mpc_clear(tau);
    for (i=0; i<6; i++) mpc_clear(thetas[i]);
    for (i=0; i<3; i++) mpc_clear(thetaconstants[i]);
    for (i=0; i<2; i++) mpc_clear(constants2[i]);

    mpfr_free_cache();

    return 0;
}
