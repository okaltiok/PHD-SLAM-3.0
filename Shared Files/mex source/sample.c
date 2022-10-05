//     This function samples from the Gaussian mixture sampling distribution
// 
//     Input:
//        obj        - struct that represent particle n of the PHD-SLAM density
//        wp         - (J x 1) GMM weights
//        mp         - (3 x J) GMM means
//        Pp         - (3 x 3 x J) GMM covariances
//        r_uniform  - uniform random variable
//        r_gaussian - 3 x 1 normal random variable
//     
//     Output:    
//        obj        - struct with updated vehicle mean and particle weight  
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

//#include "models.h"
#include "linearAlgebra.h"
#include "mex.h"
#include "matrix.h"
#include "string.h"
#include "math.h"

void loglikelihood(double *w, double *x, double *P, double *xj, int dim);

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*
     * nrhs     Number of input (right-side) arguments, or the size of the prhs array.
     * prhs     Array of input arguments.
     */
      

    /* duplicate input to output */
    plhs[0] = mxDuplicateArray(prhs[0]);    
    
    /* variable declarations of output variables */
    double* xn = mxGetDoubles(mxGetField(plhs[0], 0, "xn"));
    double* Pn = mxGetDoubles(mxGetField(plhs[0], 0, "Pn"));
    double* w = mxGetDoubles(mxGetField(plhs[0], 0, "w"));
    double* iparent = mxGetDoubles(mxGetField(plhs[0], 0, "iparent"));

    /* variable declarations of input variables */
    double* wp = mxGetDoubles(prhs[1]);      // 1 x J weight matrix
    double* mp = mxGetDoubles(prhs[2]);      // 3 x J state matrix
    double* Pp = mxGetDoubles(prhs[3]);      // 3 x 3 x J covariance matrix
    double r_uniform = mxGetScalar(prhs[4]);      // uniform random variable
    double* r_gaussian = mxGetDoubles(prhs[5]);      // 3 x 1 normal random variable
    
    /* dimension declarations */
    int J = (int)mxGetN(prhs[2]);       // number of GM components
    int xn_dim = (int)mxGetM(prhs[2]);
    int j_idx = -1;
     
    if (J > 1)
    {    
        // normalize log-weights
        double max_log_w = log(0);
        for (int j=0; j < J; j++)
        {
            if (wp[j] > max_log_w)
                max_log_w = wp[j];
        }
        
        double log_sum_w = 0.0;
        for (int j=0; j < J; j++)
        {
            log_sum_w += exp(wp[j] - max_log_w);
        }
        log_sum_w = log(log_sum_w) + max_log_w;
        
        // sample categorial distribution to obtain index j_idx
        double WP = 0.0;
        for (int j=0; j < J; j++)
        {
            wp[j] = wp[j]  - log_sum_w;
            WP += exp(wp[j]);
            if ((j_idx == -1) && (WP >= r_uniform))
            {
                j_idx = j;
            }
        }
    }
    else
    {
        j_idx = 0;
        *wp = 0.0;
    }

    // check that matrix is non-zero
    double sqrtP[9];
    if(choleskyFactorization(sqrtP, Pp + j_idx*xn_dim*xn_dim, xn_dim))
    { 
        // make a copy of the prior
        double x0[3];
        memcpy(x0, xn, xn_dim*sizeof(double));
        
        // Sample
        memcpy(xn, mp + j_idx*xn_dim, xn_dim*sizeof(double));
        matrixMultiply(xn_dim, xn_dim, xn_dim, 1, 0, 1, sqrtP, r_gaussian, xn);
        
        // add prior to particle weight
        loglikelihood(w, x0, Pn, xn, xn_dim);
        
        // compute prop and add to weight
        double prop = 0.0, L[1];
        for (int j=0; j < J; j++)
        {
            L[0] = *wp;
            loglikelihood(L, mp, Pp, xn, xn_dim);
            prop += exp(L[0]);
            
            wp ++;
            mp += xn_dim;
            Pp += xn_dim*xn_dim;
        }
        // subtract proposal from particle weight
        *w -= log(prop);
    }
            
    // set parent to zero and set elements of Pn to zero
    *iparent = 0;
    for (int j=0; j < xn_dim*xn_dim; j++)
        Pn[j] = 0.0;
}

void loglikelihood(double *w, double *x, double *P, double *xj, int dim)
{
    double Pinv[9], Pdet[1], Pnu[3], nu[3], d = 0.0;
    
    // compute residual
    for(int j = 0; j<dim; j++)
        nu[j] = xj[j] - x[j];

    // Compute d = nu'/P*nu 
    matrixInverse(Pinv, Pdet, P, dim);
    matrixMultiply(dim, dim, dim, 1, 0, 0, Pinv, nu, Pnu);
    for(int j = 0; j<dim; j++)
        d += nu[j]*Pnu[j];
    
    // compute log-likelihood
    *w -= 0.5*dim*log(2*M_PI);
    *w -= 0.5*d;
    *w -= 0.5*log(*Pdet);
}