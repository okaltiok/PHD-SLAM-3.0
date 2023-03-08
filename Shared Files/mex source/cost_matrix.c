//     This function computes the cost matrix and detection probability
// 
//     Input:
//        obj - struct that represent particle n of the PHD-SLAM density
//        Y          - {2 x m_k} matrix containing the measurements
//        params     - simulation parameters
//     
//     Output:
//        P_D        - (1 x n_k) containing the detection probabilities
//        L          - (n_k x [n_k + m_k]) cost matrix
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
#include "string.h"
#include "math.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*
     * nrhs     Number of input (right-side) arguments, or the size of the prhs array.
     * prhs     Array of input arguments.
     */
      
    /* variable declarations of field pointers */
    const mxArray* xnPtr = mxGetField(prhs[0], 0, "xn");  
    const mxArray* xlPtr = mxGetField(prhs[0], 0, "xl");  
    
    /* variable declarations of input variables */
    double* xn = mxGetDoubles(xnPtr);
    double* Pn = mxGetDoubles(mxGetField(prhs[0], 0, "Pn"));
    double* xl = mxGetDoubles(xlPtr);
    double* Pl = mxGetDoubles(mxGetField(prhs[0], 0, "Pl"));
    double* eta = mxGetDoubles(mxGetField(prhs[0], 0, "eta"));
    double* Y = mxGetDoubles(prhs[1]);      
    double* R = mxGetDoubles(mxGetField(prhs[2], 0, "R"));
    
    /* variable declarations of constants */
    double DETECTION_PROBABILITY = mxGetScalar(mxGetField(prhs[2], 0, "P_D"));
    double GATING_SIZE = mxGetScalar(mxGetField(prhs[2], 0, "gating_size"));
    double LAMBDA_C = mxGetScalar(mxGetField(prhs[2], 0, "lambda_c"));
    double FOV_RANGE = mxGetScalar(mxGetField(prhs[2], 0, "fov_range"));
    double FOV_ANGLE = mxGetScalar(mxGetField(prhs[2], 0, "fov_angle"));
   
    /* dimension declarations */
    int n_k = (int) mxGetN(xlPtr);      // number of landmarks
    int m_k = (int) mxGetN(prhs[1]);    // number of measurements
    int xn_dim = (int) mxGetScalar(mxGetField(prhs[2], 0, "xn_dim"));
    int xl_dim = (int) mxGetScalar(mxGetField(prhs[2], 0, "xl_dim"));
    int h_dim = (int) mxGetScalar(mxGetField(prhs[2], 0, "h_dim"));
    int cost_dim = m_k + n_k;
    int state_dim = xn_dim + xl_dim;

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,n_k, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_k,cost_dim,mxREAL);
    
    /* variable declarations of output variables  */
    double* PD = mxGetDoubles(plhs[0]);
    double* L = mxGetDoubles(plhs[1]);
    
    // initialize cost matrix with all -inf values
    for (int i=0; i < n_k*cost_dim; i++)
        L[i] = log(0);
    
    // help variables to store temporary data
    double mu[2], Hl[4], Hn[6], HlPl[4], HnPn[6];
    double S[4], Sinv[4], Sdet[1], Snu[2];
    double d, nu[2];
    
    /* loop through landmarks */
    for (int i=0; i < n_k; i++)
    {
        /* compute detection probability */
        adaptiveDetectionProbability(PD, xl, xn, FOV_RANGE, FOV_ANGLE, DETECTION_PROBABILITY);
        
        if (*PD > 0.0)
        {
            // compute mean and Jacobians
            h_func(xn, xl, mu, Hl, Hn);
            
            // compute S = Hl*Pl*Hl' + Hn+Pn*Hn' + R;
            memcpy(S, R, h_dim*h_dim*sizeof(double) );
            matrixMultiply(h_dim, xn_dim, xn_dim, xn_dim, 0, 0, Hn, Pn, HnPn);
            matrixMultiply(h_dim, xn_dim, h_dim, xn_dim, 1, 1, HnPn, Hn, S);
            matrixMultiply(h_dim, xl_dim, xl_dim, xl_dim, 0, 0, Hl, Pl, HlPl);
            matrixMultiply(h_dim, xl_dim, h_dim, xl_dim, 1, 1, HlPl, Hl, S);
     
            // invert S and compute determinant
            matrixInverse(Sinv, Sdet, S, h_dim);
            
            // compute constant part of log-likelihood
            *Sdet = -0.5*log(*Sdet);
            *Sdet += log(*PD) + eta[i] - 0.5*h_dim*log(2*M_PI);
            *Sdet -= log(LAMBDA_C);

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
                    L[ii*n_k] = *Sdet - 0.5 * d;
            }
        }
        
        // log-likelihood of misdetection
        L[n_k*(m_k+i)] = log(1-*PD) + eta[i];

        xl += xl_dim;
        Pl += xl_dim*xl_dim;
        L++;
        PD++;
    }
}