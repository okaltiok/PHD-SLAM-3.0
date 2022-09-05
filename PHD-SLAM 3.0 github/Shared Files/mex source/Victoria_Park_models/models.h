#ifndef MODELS_H_
#define MODELS_H_


void h_func(double *x, double *lm, double *mu, double *Hl, double *Hn);

void f_func(double *x, double *u, double dt, double *mu, double *Fx, double *Fu);

void matrixMultiply(int rA, int cA, int rB, int cB, int transpose, int opt, double *A, double *B, double *C);

void matrixInverse(double *invX, double *det, double *X, int dim);

void choleskyDecomposition(double *P, double *L);

double pi_to_pi(double angle);

void adaptiveDetectionProbability(double *PD, double *xl, double *xn, double FOV_RANGE, double FOV_ANGLE, double DETECTION_PROBABILITY);


#endif