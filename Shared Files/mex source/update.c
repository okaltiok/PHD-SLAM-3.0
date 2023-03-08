//     This function computes the update step of PHD-SLAM  filter
// 
//     Input:
//        obj        - struct that represent particle n of the PHD-SLAM density at time k | k-1
//        y          - {2 x m_k} matrix containing the measurements
//        params     - simulation parameters
//     
//     Output:
//        obj        - struct that represent particle n of the PHD-SLAM density at time k | k
//     
//     Author   : Ossi Kaltiokallio
//                Tampere University, Department of Electronics and
//                Communications Engineering
//                Korkeakoulunkatu 1, 33720 Tampere
//                ossi.kaltiokallio@tuni.fi
//     Last Rev : 26/8/2022
//     Tested   : '9.8.0.1359463 (R2020a) Update 1'
//     
//     Copyright notice: You are free to modify, extend and distribute
//        this code granted that the author of the original code is
//        mentioned as the original author of the code.


#include "models.h"
#include "linearAlgebra.h"
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include "string.h"
#include "math.h"



/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*
     * nrhs     Number of input (right-side) arguments, or the size of the prhs array.
     * prhs     Array of input arguments.
     */
      
    /* duplicate input to output */
    plhs[0] = mxDuplicateArray(prhs[0]);
        
    int m_k = (int) mxGetN(prhs[1]);      // number of measurements
    if(m_k == 0)
    {
        return;
    }
    
    /* variable declarations of field pointers */
    mxArray* xl_ptr = mxGetField(plhs[0], 0, "xl");
    mxArray* Pl_ptr = mxGetField(plhs[0], 0, "Pl"); 
    mxArray* eta_ptr = mxGetField(plhs[0], 0, "eta"); 
    mxArray* w_ptr = mxGetField(plhs[0], 0, "w"); 
    mxArray* birth_ptr = mxGetField(plhs[0], 0, "birth_y"); 
         
    /* variable declarations of output variables */
    double* xl = mxGetDoubles(xl_ptr);
    double* Pl = mxGetDoubles(Pl_ptr);
    double* eta = mxGetDoubles(eta_ptr);
    double* w = mxGetDoubles(w_ptr);
    double* birth_y = mxGetDoubles(birth_ptr);
    
    /* variable declarations of input variables */
    double* xn = mxGetDoubles(mxGetField(prhs[0], 0, "xn"));
    double* Y = mxGetDoubles(prhs[1]);      
    double* R = mxGetDoubles(mxGetField(prhs[2], 0, "R"));
    
    /* variable declarations of constants */
    double DETECTION_PROBABILITY = mxGetScalar(mxGetField(prhs[2], 0, "P_D"));
    double GATING_SIZE = mxGetScalar(mxGetField(prhs[2], 0, "gating_size"));
    double LAMBDA_C = mxGetScalar(mxGetField(prhs[2], 0, "lambda_c"));
    double FOV_RANGE = mxGetScalar(mxGetField(prhs[2], 0, "fov_range"));
    double FOV_ANGLE = mxGetScalar(mxGetField(prhs[2], 0, "fov_angle"));
    
    /* dimension declarations */
    int n_k = (int) mxGetN(xl_ptr);      // number of landmarks
    int xn_dim = (int) mxGetScalar(mxGetField(prhs[2], 0, "xn_dim"));
    int xl_dim = (int) mxGetScalar(mxGetField(prhs[2], 0, "xl_dim"));
    int h_dim = (int) mxGetScalar(mxGetField(prhs[2], 0, "h_dim"));

    
    /* allocate memory for temporary variables used in the update */
    double* PD = mxCalloc(n_k*sizeof(double),sizeof(double));
    double* hh = mxCalloc(h_dim*n_k*sizeof(double),sizeof(double));
    double* L = mxCalloc(n_k*m_k*sizeof(double),sizeof(double));
    double* KK = mxCalloc(xl_dim*h_dim*n_k*sizeof(double),sizeof(double));
    double* PP = mxCalloc(xl_dim*xl_dim*n_k*sizeof(double),sizeof(double));
    double* eta_tilde_sum = mxCalloc(m_k*sizeof(double),sizeof(double));
    int* DA = mxCalloc(m_k*sizeof(int),sizeof(int));
    int NDA = 0;
    
    // initialize cost matrix with all -inf values
    for (int i=0; i < n_k*m_k; i++)
        L[i] = log(0);

    // help variables
    double mu[2], Hl[4], Hn[6], HlPl[4], PlHl[4];
    double S[4], Sinv[4], Sdet[1], Snu[2], KS[4];
    double d, nu[2];
    
    /* loop through landmarks and measurements and compute eta, Kalman gain 
     * and covariance used in the PHD update */
    for (int i=0; i < n_k; i++)
    {
        /* compute detection probability */
        adaptiveDetectionProbability(PD, xl, xn, FOV_RANGE, FOV_ANGLE, DETECTION_PROBABILITY);
        
        if (*PD > 0.0)
        {
            memcpy(PP, Pl, xl_dim*xl_dim*sizeof(double));

            // compute mean and Jacobians
            h_func(xn, xl, mu, Hl, Hn);

            // compute S = Hl*Pl*Hl' + R;
            memcpy(S, R, h_dim*h_dim*sizeof(double) );
            matrixMultiply(h_dim, xl_dim, xl_dim, xl_dim, 0, 0, Hl, Pl, HlPl);
            matrixMultiply(h_dim, xl_dim, h_dim, xl_dim, 1, 1, HlPl, Hl, S);
            
            // invert S and compute determinant
            matrixInverse(Sinv, Sdet, S, h_dim);
                        
            // compute constant part of log-likelihood
            *Sdet = -0.5*log(*Sdet);
            *Sdet += log(*PD) + eta[i] - 0.5*h_dim*log(2*M_PI);

            // compute log-likelihoods
            for(int ii = 0; ii<m_k; ii++)
            {
                // compute residual
                nu[0] = Y[ii*h_dim] - mu[0];
                nu[1] = pi_to_pi(Y[ii*h_dim+1] - mu[1]);

                /* compute d = nu'/S*nu */
                matrixMultiply(h_dim, h_dim, h_dim, 1, 0, 0, Sinv, nu, Snu);
                d = nu[0]*Snu[0] + nu[1]*Snu[1];
                
                // perform gating and store log-likelihood of detection
                if (d < GATING_SIZE)
                {
                    DA[ii] += 1;
                    NDA++;
                    L[ii*n_k] = *Sdet - 0.5 * d;
                    eta_tilde_sum[ii] += exp(*(L + ii*n_k));
                }
            }

             /* compute K = Pl*Hl'*Sinv */
            matrixMultiply(xl_dim, xl_dim, h_dim, xl_dim, 1, 0, PP, Hl, PlHl);
            matrixMultiply(h_dim, xl_dim, h_dim, h_dim, 0, 0, PlHl, Sinv, KK);
            
            /* compute Pl = Pl - K*S*K' */
            matrixMultiply(xl_dim, h_dim, h_dim, h_dim, 0, 0, KK, S, KS);
            matrixMultiply(xl_dim, h_dim, xl_dim, h_dim, 1, -1, KS, KK, PP);

            // copy output variables
            memcpy(hh, mu,h_dim*sizeof(double));
        }
        hh += h_dim;
        KK += h_dim*xl_dim;
        PP += xl_dim*xl_dim;
        xl += xl_dim;
        Pl += xl_dim*xl_dim;
        PD++;
        L++;
    }
    
    /* move pointers to beginning */
    hh -= h_dim*n_k;
    KK -= h_dim*xl_dim*n_k;
    PP -= xl_dim*xl_dim*n_k;
    xl -= xl_dim*n_k;
    Pl -= xl_dim*xl_dim*n_k;
    PD -= n_k;
    L -= n_k;

    /* determine number of birth components */
    int n_birth = 0;
    
    for(int ii = 0; ii<m_k; ii++)
    {
        if (DA[ii] == 0)
            n_birth++;
    }
    
    /* allocate memory for birth mesurements */
    if (n_birth > 0)
    {
        const mwSize birth_ndim = 2;
        const mwSize birth_dims[2] = {h_dim, n_birth};
        mxSetDimensions(birth_ptr, birth_dims, birth_ndim);
        birth_y = mxRealloc(birth_y,n_birth*h_dim*sizeof(double));
        mxSetPr(birth_ptr,birth_y);
    }
    
    /* define dimensions of landmark parameters */
    mwSize xl_dims[2] = {xl_dim,n_k + NDA};
    mwSize Pl_dims[3] = {xl_dim,xl_dim,n_k + NDA};
    mwSize eta_dims[2] = {1,n_k + NDA};
        
    if (n_k == 0)
    {
        xl_dims[0] = 0;
        Pl_dims[0] = 0;
        Pl_dims[1] = 0;
        eta_dims[0] = 0;
    }
    
    /* reallocate memory for landmarks */
    mxSetDimensions(xl_ptr, xl_dims, 2);
    xl = mxRealloc(xl,(n_k+NDA)*xl_dim*sizeof(double));
    mxSetPr(xl_ptr,xl);
    
    mxSetDimensions(eta_ptr, eta_dims, 2);
    eta = mxRealloc(eta,(n_k+NDA)*sizeof(double));
    mxSetPr(eta_ptr,eta);

    mxSetDimensions(Pl_ptr, Pl_dims, 3);
    Pl = mxRealloc(Pl,(n_k+NDA)*xl_dim*xl_dim*sizeof(double));
    mxSetPr(Pl_ptr,Pl);


    /* update parameters of misdetection and move pointers*/
    for (int i=0; i < n_k; i++)
    {
        eta[i] = log(1 - PD[i]) + eta[i];
    }
    eta += n_k;
    Pl += xl_dim*xl_dim*n_k;

    double tmp;
    int birth_flag;
    int k = n_k;
    
    // update particle weight and parameters of detected landmarks
    for (int i=0; i < m_k; i++)
    {
        tmp = log(LAMBDA_C + eta_tilde_sum[i]);
        *w += tmp;
        birth_flag = 1;
        
        for (int j=0; j < n_k; j++)
        {
            if (L[i*n_k + j] > log(0))
            {
                // compute residual
                nu[0] = Y[i*h_dim] - hh[j*h_dim];
                nu[1] = pi_to_pi(Y[i*h_dim+1] - hh[j*h_dim+1]);
                
                /* compute xl = xl + K*nu */
                memcpy(xl+xl_dim*k, xl+xl_dim*j, xl_dim*sizeof(double) );
                matrixMultiply(xl_dim, h_dim, h_dim, 1, 0, 1, KK+xl_dim*h_dim*j, nu, xl+xl_dim*k);
                
                /* update Pl, eta */
                memcpy(Pl, PP+xl_dim*xl_dim*j, xl_dim*xl_dim*sizeof(double) );
                *eta = L[i*n_k + j] - tmp;
                
                k++;
                eta++;
                Pl += xl_dim*xl_dim;
                birth_flag = 0;
            }
        }
        
        /* 
         * measurement isn't associated to any landmark, append measurement
         * to create a birth during the PHD prediction step
         */
        if (birth_flag == 1)
        {
            memcpy(birth_y, Y + i*h_dim, h_dim*sizeof(double) );
            birth_y += h_dim;
        }
    }
    
    // free allocated memory
    mxFree(PD);
    mxFree(hh);
    mxFree(L);
    mxFree(KK);
    mxFree(PP);
    mxFree(eta_tilde_sum);
    mxFree(DA);
}