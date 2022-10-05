//     This function computes the partitioned IPL-OID approximation for J
//     Gaussian mixture components
// 
//     Input:
//        obj     - struct that represent particle n of the PHD-SLAM density
//        col4row - (n_k x J) assignment matrix 
//        Y       - {2 x m_k} measurement matrix
//        P_D     - (1 x n_k) detection probabilities
//        params  - simulation parameters
//     
//     Output:
//        WPI      - (J x 1) GMM weights
//        MPI      - (3 x J) GMM means
//        PPI      - (3 x 3 x J) GMM means
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

static const ptrdiff_t mm = 2;
static const ptrdiff_t nn = 3;
static const ptrdiff_t nl = 2;
   
#define XN_DIM_SOD 24
#define XN_DIM2_SOD 72
#define H_DIM_SOD 16
#define H_DIM2_SOD 32
#define XNH_DIM_SOD 48

/* function declarations */
void moments(double *xn, double *Pn, double *xl, double *Pl, double *mu, double *Pys, double *S);

void linearization(double *mpi, double *Ppi, double *invPpi, double *Pys, double *A, double *b, double *Omega);

void update(double *ms, double *Ps, double *A, double *b,  double *Omega,  double *y, double *gamma);

void KL_divergence(double *dkls, double *Ps, double *invPs, double *detPs, double *ms, double *Ppi, double *detPpi, double *mpi);

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
    double* msp = mxGetDoubles(xnPtr);
    double* Psp = mxGetDoubles(mxGetField(prhs[0], 0, "Pn"));
    double* xl = mxGetDoubles(xlPtr);
    double* Pl = mxGetDoubles(mxGetField(prhs[0], 0, "Pl"));
    double* eta = mxGetDoubles(mxGetField(prhs[0], 0, "eta"));
    double* col4row = mxGetDoubles(prhs[1]);      // assignment matrix   
    double* Y = mxGetDoubles(prhs[2]);      // 2 x m_k measurement matrix
    double* P_D = mxGetDoubles(prhs[3]);      // 1 x n_k detection probabilities
    double* R = mxGetDoubles(mxGetField(prhs[4], 0, "R"));
    
    /* variable declarations of constants */
    int L = mxGetScalar(mxGetField(prhs[4], 0, "L"));
    double epsilon = mxGetScalar(mxGetField(prhs[4], 0, "epsilon"));
    double* gamma_table = mxGetDoubles(mxGetField(prhs[4], 0, "gamma"));
    double lambda_c = mxGetScalar(mxGetField(prhs[4], 0, "lambda_c"));

    /* dimension declarations */
    int J = (int)mxGetN(prhs[1]);       // number of GM components
    int n_k = (int) mxGetN(xlPtr);      // number of landmarks
    int m_k = (int) mxGetN(prhs[2]);    // number of landmarks
    int xn_dim = (int) mxGetScalar(mxGetField(prhs[4], 0, "xn_dim"));
    int xl_dim = (int) mxGetScalar(mxGetField(prhs[4], 0, "xl_dim"));
    int h_dim = (int) mxGetScalar(mxGetField(prhs[4], 0, "h_dim"));
    int xn_dim2 = xn_dim*xn_dim;
    int h_dim2 = h_dim*h_dim;
    
    /* create the output matrix */
    const mwSize ndim = 3;
    const mwSize dims3D[3] = {xn_dim,xn_dim,J};
    
    plhs[0] = mxCreateDoubleMatrix(1, J, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(xn_dim, J, mxREAL);
    plhs[2] = mxCreateNumericArray(ndim, dims3D, mxDOUBLE_CLASS, mxREAL);
    
    double* WPI = mxGetDoubles(plhs[0]);
    double* MPI = mxGetDoubles(plhs[1]);
    double* PPI = mxGetDoubles(plhs[2]);

    // help variables
    double wpi[1], mpi[3], Ppi[9];
    double ms[3], Ps[9], y[2], my[2], Py[4], Pys[6];
    double invPpi[9], detPpi[1], invPs[9], detPs[1];
    double A[6], b[2], Omega[4];
    double nu[2], S[4], invS[4], Snu[2];
    double gamma[1], dkls[1];

    int done, l;
    double gammaT;

    // compute IPL-OID approximation for each GMM component
    for (int j=0; j < J; j++)
    {
        // obtain gamma threshold from a lookup table
        l = 0;
        for (int i=0; i < n_k; i++)
        {
            if ((int)col4row[i] <= m_k)
            {
                l++;
            }
        }
        
        gammaT = gamma_table[l];

        // initialize linearization density
        memcpy(mpi, msp, XN_DIM_SOD);
        memcpy(Ppi, Psp, XN_DIM2_SOD);

        done = 0;
        l = 0;

        while (done == 0)
        {
            l++;
            *gamma = 0.0;
            
            /* copy priors to ms and Ps */
            memcpy(ms, msp, XN_DIM_SOD);
            memcpy(Ps, Psp, XN_DIM2_SOD);
            
            // compute inverse of linearization covariance
            matrixInverse(invPpi, detPpi, Ppi, xn_dim);
            
            /* compute partitioned OID approximation for l-th iteration */
            for (int i=0; i < n_k; i++)
            {
                if ((int)col4row[i] <= m_k)
                {
                    /* Compute moments */
                    memcpy(Omega, R, H_DIM2_SOD);
                    moments(mpi, Ppi, xl, Pl, b, Pys, Omega);
                    
                    /* Compute linearization */
                    linearization(mpi, Ppi, invPpi, Pys, A, b, Omega);
                    
                    /*update ms, Ps and compute test statistic gamma */
                    update(ms, Ps, A, b, Omega,  Y + (int)((col4row[i]-1)*h_dim), gamma);
                }   
                xl += xl_dim;
                Pl += xl_dim*xl_dim;
            }
            

            // check if posterior update was successful, if not, exit loop
            // and use the previous best approximation
            if (*gamma >= gammaT)
            {
                done = 1;
            }
            else
            {
                if(l > 1)
                {
                    *dkls = 0.0;
                    
                    // compute inverse
                    matrixInverse(invPs, detPs, Ps, xn_dim);

                    // compute change in KL divergence
                    KL_divergence(dkls,Ps,invPs,detPs,ms,Ppi,detPpi,mpi);
                    if(*dkls < epsilon)
                    {
                        done = 1;
                    }
                }

                memcpy(mpi, ms, XN_DIM_SOD);
                memcpy(Ppi, Ps, XN_DIM2_SOD);
            }
            
            // Convergence criteria and tests
            if (l >= L)
            {
                done = 1;
            }

            xl -= xl_dim*n_k;
            Pl -= xl_dim*xl_dim*n_k;
        }
    
        /* IPL-OID approximation complete, compute weight for j-th component */
        *wpi = 0.0;
        if(J > 1)
        {
            for (int i=0; i < n_k; i++)
            {
                if ((int)col4row[i] <= m_k)
                {    
                    memcpy(y, Y + (int)((col4row[i]-1)*h_dim), H_DIM_SOD );
                    memcpy(Py, R, H_DIM2_SOD );
                    
                    // Compute moments 
                    moments(mpi, Ppi, xl, Pl, my, Pys, Py);
                    
                    // invert Py and compute determinant
                    matrixInverse(invS, detPs, Py, h_dim);
                    
                    // compute residual
                    nu[0] = y[0] - my[0];
                    nu[1] = pi_to_pi(y[1] - my[1]);
                    
                    /* compute gamma = nu'/S*nu */
                    matrixMultiply(mm, mm, mm, 1, 0, 0, invS, nu, Snu);
                    *gamma = nu[0]*Snu[0] + nu[1]*Snu[1];

                    // update log-weight with detection
                    *wpi += P_D[i] + eta[i] - log(lambda_c);
                    *wpi -= 0.5*(h_dim*log(2*M_PI) + log(*detPs) + *gamma);
                }
                else
                {
                    // update log-weight with misdetection
                    *wpi += log(1 - exp(P_D[i])) + eta[i];
                }
                xl += xl_dim;
                Pl += xl_dim*xl_dim;
            }
            xl -= xl_dim*n_k;
            Pl -= xl_dim*xl_dim*n_k;
        }

        // copy GMM parameters to output 
        memcpy(WPI, wpi, sizeof(double));
        memcpy(MPI, mpi, XN_DIM_SOD);
        memcpy(PPI, Ppi, XN_DIM2_SOD);

        WPI++;
        MPI += xn_dim;
        PPI += xn_dim2;
        col4row += n_k;
    }
}

void KL_divergence(double *dkls, double *Ps, double *invPs, double *detPs, double *ms, double *Ppi, double *detPpi, double *mpi)
{   
    double nu[3], traceP[9], Pnu[3]; 
        
    /* compute trace(Ps\Ppi)  */
    matrixMultiply(nn, nn, nn, nn, 0, 0, invPs, Ppi, traceP);
    
    for(int j = 0; j < nn; j++)
    {
        nu[j] = ms[j] - mpi[j];
        *dkls += *(traceP + j*nn+j);
    }
    
    *dkls -= log(*detPpi);
    *dkls += log(*detPs);
    *dkls -= nn;
    
    matrixMultiply(nn, nn, nn, 1, 0, 0, invPs, nu, Pnu);
    for(int j = 0; j < nn; j++)
    {
        *dkls += nu[j]*Pnu[j];
    }
    
    *dkls /= 2.0;
}

void update(double *ms, double *Ps, double *A, double *b,  double *Omega,  double *y, double *gamma)
{    
    /* 
     * Update:
     * S = A*Ps*A' + Omega
     * K = Ps*A'/S
     * ms = ms + K*nu
     * Ps = Ps - K*S*K' 
     */
    
    double AP[6], PsA[6], K[6], KS[6], invS[4], detS[1], nu[2], Snu[2];
    
    /* compute S = A*Ps*A' + Omega */
    matrixMultiply(mm, nn, nn, nn, 0, 0, A, Ps, AP);
    matrixMultiply(mm, nn, mm, nn, 1, 1, AP, A, Omega);
    
    // invert matrix
    matrixInverse(invS, detS, Omega, mm);

    /* compute K = Ps*A'*invS */
    matrixMultiply(nn, nn, mm, nn, 1, 0, Ps, A, PsA);
    matrixMultiply(nn, mm, mm, mm, 0, 0, PsA, invS, K);
    
    /* compute b = A*ms + b* */
    matrixMultiply(mm, nn, nn, 1, 0, 1, A, ms, b);
    
    /* compute nu = y - A*ms - b* */
    nu[0] = y[0] - b[0];
    nu[1] = pi_to_pi(y[1] - b[1]);
    
    /* compute ms = ms + K*nu */
    matrixMultiply(nn, mm, mm, 1, 0, 1, K, nu, ms);
    
    /* compute P = P - K*S*K' */
    matrixMultiply(nn, mm, mm, mm, 0, 0, K, Omega, KS);
    matrixMultiply(nn, mm, nn, mm, 1, -1, KS, K, Ps);
    
    /* compute gamma += nu'/S*nu */
    matrixMultiply(mm, mm, mm, 1, 0, 0, invS, nu, Snu);
    *gamma += nu[0]*Snu[0] + nu[1]*Snu[1];
    
}


void linearization(double *mpi, double *Ppi, double *invPpi, double *Pys, double *A, double *b, double *Omega)
{
    /* 
     * Linearization:
     * A = Pys/Ppi
     * b = my - A*mpi
     * Omega = Py - A*Ppi*A' 
     */
    
    double AP[6];

    /* compute A = Pys/Ppi */
    matrixMultiply(mm, nn, nn, nn, 0, 0, Pys, invPpi, A);

    /* compute b = my - A*mpi */
    matrixMultiply(mm, nn, nn, 1, 0, -1, A, mpi, b);

    /* compute Omega = Py - A*Ppi*A' */
    matrixMultiply(mm, nn, nn, nn, 0, 0, A, Ppi, AP);
    matrixMultiply(mm, nn, mm, nn, 1, -1, AP, A, Omega);
}        

void moments(double *xn, double *Pn, double *xl, double *Pl, double *mu, double *Pys, double *Py)
{
    /* 
     *  Moments:
     *  mu = g(xn,xl)
     *  Py = Hn*Pn*Hn' + Hl*Pl*Hl' + R
     */
    
    double Hl[4], Hn[6], HlPl[4];
            
    /* compute mean and Jacobians */
    h_func(xn, xl, mu, Hl, Hn);
   
    /* 1st, compute Py += Hn*Pn*Hn' */
    matrixMultiply(mm, nn, nn, nn, 0, 0, Hn, Pn, Pys);
    matrixMultiply(mm, nn, mm, nn, 1, 1, Pys, Hn, Py);
    
    /* 2nd, compute Py += Hl*Pl*Hl'*/
    matrixMultiply(mm, nl, nl, nl, 0, 0, Hl, Pl, HlPl);
    matrixMultiply(mm, nl, mm, nl, 1, 1, HlPl, Hl, Py);
}
                    
