#ifndef MODELS_H_
#define MODELS_H_

void h_func(double *x, double *lm, double *mu, double *Hl, double *Hn);

void f_func(double *x, double *u, double dt, double *mu, double *Fx, double *Fu);

void adaptiveDetectionProbability(double *PD, double *xl, double *xn, double FOV_RANGE, double FOV_ANGLE, double DETECTION_PROBABILITY);

#endif