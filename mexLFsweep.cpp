#include <math.h>
#include "matrix.h"
#include "mex.h"   //--This one is required
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

double sweepingLF(double**** phi, int N1, int N2, int N3, int N4, double***** xs, double* dx, double a_upper, double a_lower, double omega_max, double**** obs)
{
    int i,j,k,h;
    int s1,s2,s3,s4;
    double p, q, r, t;
    double c, H, phiTemp, phiOld, error;
    double obs_this;
    double sigma1, sigma2, sigma3, sigma4;
    double v;
    double a, omega; // optimal control (a, omega)

    
    double tmp;
    
    int count;
    count = 0;
//     double small=0.1;
//     double pP, pM, qP, qM, rP, rM, HP, HM, sigmaP, sigmaM;
    
    error = 0;
//     N1 = N[0]; N2 = N[1]; N3 = N[2];
    
    // mexPrintf("here in sweeping beginning place");
    for (s1=1; s1>=-1; s1-=2 )
    for (s2=1; s2>=-1; s2-=2 )
    for (s3=1; s3>=-1; s3-=2 )
    for (s4=1; s4>=-1; s4-=2 )
    {
        // LF sweeping module
        for ( i=(s1<0 ? (N1-2):1); (s1<0 ? i>=1:i<=(N1-2)); i+=s1 )
        for ( j=(s2<0 ? (N2-2):1); (s2<0 ? j>=1:j<=(N2-2)); j+=s2 ) 
        for ( k=(s3<0 ? (N3-1):0); (s3<0 ? k>=0:k<=(N3-1)); k+=s3 ) 
        for ( h=(s4<0 ? (N4-2):1); (s4<0 ? h>=1:h<=(N4-2)); h+=s4 )
        {
            phiOld = phi[i][j][k][h];
            obs_this = obs[i][j][k][h];
            // mexPrintf("phiold in front:%f", phiOld);
            
            count += 1;
//             c = (sigma1[i][j][k]/dx[0] + sigma2[i][j][k]/dx[1] + sigma3/dx[2]);
            
            if (k == 0)
            {
                p = (phi[i+1][j][k][h] - phi[i-1][j][k][h])/(2*dx[0]);  // dVx
                q = (phi[i][j+1][k][h] - phi[i][j-1][k][h])/(2*dx[1]);  // dVy
                r = (phi[i][j][k+1][h] - phi[i][j][N3-1][h])/(2*dx[2]); // dVtheta
                t = (phi[i][j][k][h+1] - phi[i][j][k][h-1])/(2*dx[3]);  // dVv
//                 
//                 if (obs_this > 0) {
//                     if (t <= 0) { //TODO?
//                         a = 0;
//                     }
//                     // v = 0;
//                     // omega = 0;
//                     // d = 0;
//                 } else {
                    // Compute optimal control
                    if (t >= 0) { //TODO?
                        a = a_lower;
                    } else { 
                        a = a_upper;
                    } 
                    if (r >= 0) { //TODO?
                        omega = -omega_max;
                    } else { 
                        omega = omega_max;
                    }
//                 }
//                 sigma1 =  abs(-transDrag*xs[i][j][k][h][0]/m + T1*sin(xs[i][j][k][h][2])/m + T2*sin(xs[i][j][k][h][2])/m);
//                 sigma2 =  abs(-grav - transDrag*xs[i][j][k][h][1]/m + T1*cos(xs[i][j][k][h][2])/m + T2*cos(xs[i][j][k][h][2])/m);
//                 sigma3 =  abs(xs[i][j][k][h][3]);
//                 sigma4 =  abs(-rotDrag*xs[i][j][k][h][3]/Iyy + L*T1/Iyy - L*T2/Iyy);
                
//                 v = sqrt(pow(xs[i][j][k][h][1],2)+pow(xs[i][j][k][h][2],2));
                v = xs[i][j][k][h][3]; //is this correct?
//                 v = 1;
                sigma1 = abs(v*cos(xs[i][j][k][h][2])); // x1_dot = x4*cos(x3)
                sigma2 = abs(v*sin(xs[i][j][k][h][2])); // x2_dot = x4*sin(x3)
                sigma3 = abs(omega);                    // x3_dot = u1
                sigma4 = abs(a);                        // x4_dot = u2

                c = (sigma1/dx[0] + sigma2/dx[1] + sigma3/dx[2] + sigma4/dx[3]);
                
                // H = p*F
                H = (-1) * (v*cos(xs[i][j][k][h][2])*p + v*sin(xs[i][j][k][h][2])*q + omega*r +a*t +1); //TODO? in4d still correct?
           
                phiTemp = -H + sigma1*(phi[i+1][j][k][h] + phi[i-1][j][k][h])/(2*dx[0])\
                        + sigma2*(phi[i][j+1][k][h] + phi[i][j-1][k][h])/(2*dx[1])\
                        + sigma3*(phi[i][j][k+1][h] + phi[i][j][N3-1][h])/(2*dx[2])\
                        + sigma4*(phi[i][j][k][h+1] + phi[i][j][k][h-1])/(2*dx[3]);

                }
            else if (k == N3-1)
            {
//                 p = (phi[i+1][j][k] - phi[i-1][j][k])/(2*dx[0]);
//                 q = (phi[i][j+1][k] - phi[i][j-1][k])/(2*dx[1]);
//                 r = (phi[i][j][0] - phi[i][j][k-1])/(2*dx[2]);
                
                p = (phi[i+1][j][k][h] - phi[i-1][j][k][h])/(2*dx[0]);
                q = (phi[i][j+1][k][h] - phi[i][j-1][k][h])/(2*dx[1]);
                r = (phi[i][j][0][h] - phi[i][j][k-1][h])/(2*dx[2]);
                t = (phi[i][j][k][h+1] - phi[i][j][k][h-1])/(2*dx[3]);
                
//                 if (obs_this > 0) {
//                     if (t <= 0) { //TODO?
//                         a = 0;
//                     }
//                 } else {
                    // Compute optimal control
                     if (t >= 0) { //TODO?
                        a = a_lower;
                    } else { 
                        a = a_upper;
                    }
                    if (r >= 0) { //TODO?
                        omega = - omega_max;
                    } else { 
                        omega = omega_max;
                    }
//                 }
                
//                 sigma1 =  abs(-transDrag*xs[i][j][k][h][0]/m + T1*sin(xs[i][j][k][h][2])/m + T2*sin(xs[i][j][k][h][2])/m);
//                 sigma2 =  abs(-grav - transDrag*xs[i][j][k][h][1]/m + T1*cos(xs[i][j][k][h][2])/m + T2*cos(xs[i][j][k][h][2])/m);
//                 sigma3 =  abs(xs[i][j][k][h][3]);
//                 sigma4 =  abs(-rotDrag*xs[i][j][k][h][3]/Iyy + L*T1/Iyy - L*T2/Iyy);
                
//                 v = sqrt(pow(xs[i][j][k][h][1],2)+pow(xs[i][j][k][h][2],2));
                v = xs[i][j][k][h][3];
//                 v = 1;
                sigma1 = abs(v*cos(xs[i][j][k][h][2])); // x1_dot = x4*cos(x3)
                sigma2 = abs(v*sin(xs[i][j][k][h][2])); // x2_dot = x4*sin(x3)
                sigma3 = abs(omega);                    // x3_dot = u1
                sigma4 = abs(a);                        // x4_dot = u2

                c = (sigma1/dx[0] + sigma2/dx[1] + sigma3/dx[2] + sigma4/dx[3]);
                
                // H = p*F
                H = (-1) * (v*cos(xs[i][j][k][h][2])*p + v*sin(xs[i][j][k][h][2])*q + omega*r +a*t +1); //TODO?
                
                phiTemp = -H + sigma1*(phi[i+1][j][k][h] + phi[i-1][j][k][h])/(2*dx[0])\
                        + sigma2*(phi[i][j+1][k][h] + phi[i][j-1][k][h])/(2*dx[1])\
                        + sigma3*(phi[i][j][0][h] + phi[i][j][k-1][h])/(2*dx[2])\
                        + sigma4*(phi[i][j][k][h+1] + phi[i][j][k][h-1])/(2*dx[3]);
            }
            else
            {
//                 p = (phi[i+1][j][k] - phi[i-1][j][k])/(2*dx[0]);
//                 q = (phi[i][j+1][k] - phi[i][j-1][k])/(2*dx[1]);
//                 r = (phi[i][j][k+1] - phi[i][j][k-1])/(2*dx[2]);
                
                p = (phi[i+1][j][k][h] - phi[i-1][j][k][h])/(2*dx[0]);
                q = (phi[i][j+1][k][h] - phi[i][j-1][k][h])/(2*dx[1]);
                r = (phi[i][j][k+1][h] - phi[i][j][k-1][h])/(2*dx[2]);
                t = (phi[i][j][k][h+1] - phi[i][j][k][h-1])/(2*dx[3]);
                
//                 if (obs_this > 0) {
//                     if (t <= 0) { //TODO?
//                         a = 0;
//                     }
//                 } else {
                // Compute optimal control
                     if (t >= 0) { //TODO?
                        a = a_lower;
                    } else { 
                        a = a_upper;
                    } 
                    if (r >= 0) { //TODO?
                        omega = - omega_max;
                    } else { 
                        omega = omega_max;
                    }
//                 }
//                 
//                 sigma1 =  abs(-transDrag*xs[i][j][k][h][0]/m + T1*sin(xs[i][j][k][h][2])/m + T2*sin(xs[i][j][k][h][2])/m);
//                 sigma2 =  abs(-grav - transDrag*xs[i][j][k][h][1]/m + T1*cos(xs[i][j][k][h][2])/m + T2*cos(xs[i][j][k][h][2])/m);
//                 sigma3 =  abs(xs[i][j][k][h][3]);
//                 sigma4 =  abs(-rotDrag*xs[i][j][k][h][3]/Iyy + L*T1/Iyy - L*T2/Iyy);
                
//                 v = 1;
//                 v = sqrt(pow(xs[i][j][k][h][1],2)+pow(xs[i][j][k][h][2],2));
                v = xs[i][j][k][h][3];

                sigma1 = abs(v*cos(xs[i][j][k][h][2])); // x1_dot = x4*cos(x3)
                sigma2 = abs(v*sin(xs[i][j][k][h][2])); // x2_dot = x4*sin(x3)
                sigma3 = abs(omega);                    // x3_dot = u1
                sigma4 = abs(a);                        // x4_dot = u2

                c = (sigma1/dx[0] + sigma2/dx[1] + sigma3/dx[2] + sigma4/dx[3]);
                
                // H = p*F
                H = (-1) * (v*cos(xs[i][j][k][h][2])*p + v*sin(xs[i][j][k][h][2])*q + omega*r +a*t +1); //TODO?                
                
                // phiTemp = - H + sigma1*(phi[i+1][j][k] + phi[i-1][j][k])/(2*dx[0]) + sigma2*(phi[i][j+1][k] + phi[i][j-1][k])/(2*dx[1]) + sigma3*(phi[i][j][k+1] + phi[i][j][k-1])/(2*dx[2]);
                phiTemp = -H + sigma1*(phi[i+1][j][k][h] + phi[i-1][j][k][h])/(2*dx[0])\
                        + sigma2*(phi[i][j+1][k][h] + phi[i][j-1][k][h])/(2*dx[1])\
                        + sigma3*(phi[i][j][k+1][h] + phi[i][j][k-1][h])/(2*dx[2])\
                        + sigma4*(phi[i][j][k][h+1] + phi[i][j][k][h-1])/(2*dx[3]);
            }
            
            if (obs_this > 0) {
                phi[i][j][k][h] = 100;
//                 phi[i][j][k][h] = max(phiTemp/c, phiOld);
            }else{
                phi[i][j][k][h] = min(phiTemp/c, phiOld);
            }
            error = max(error, phiOld - phi[i][j][k][h]);
//             mexPrintf("phitemp:%10f, phiold:%f \n", (phiTemp, phiOld));
//             mexPrintf("c:%f \n", c);
//             mexPrintf("phitemp/c:%10f \n", phiTemp/c);
//             
//             tmp = min(phiTemp/c, phiOld);
//             mexPrintf("tmp:%10f", tmp);
//             phi[i][j][k][h] = tmp;
//             error = max(error, phiOld - phi[i][j][k][h]);
//             if (count >= 5)
//             {
//                 return 0;
//             }
            
            

        }
        // computational boundary condition
        // mexPrintf("here we finish one sweeping module");
        for ( j = 0; j <= (N2-1); j++)
        for ( k = 0; k <= (N3-1); k++)
        for ( h = 0; h <= (N4-1); h++)
        {
            phiOld = phi[0][j][k][h]; 
            phi[0][j][k][h] = min(max(2*phi[1][j][k][h] - phi[2][j][k][h], phi[2][j][k][h]), phiOld); 
            error = max(error, phiOld - phi[0][j][k][h]);
            
            phiOld = phi[N1-1][j][k][h];
            phi[N1-1][j][k][h] = min(max(2*phi[N1-2][j][k][h] - phi[N1-3][j][k][h], phi[N1-3][j][k][h]), phiOld); 
            error = max(error, phiOld - phi[N1-1][j][k][h]);
        }

        for ( k = 0; k <= (N3-1); k++)
        for ( i = 0; i <= (N1-1); i++)
        for ( h = 0; h <= (N4-1); h++)
        {
            phiOld = phi[i][0][k][h]; 
            phi[i][0][k][h] = min(max(2*phi[i][1][k][h] - phi[i][2][k][h], phi[i][2][k][h]), phiOld); 
            error = max(error, phiOld - phi[i][0][k][h]);
            
            phiOld = phi[i][N2-1][k][h];
            phi[i][N2-1][k][h] = min(max(2*phi[i][N2-2][k][h] - phi[i][N2-3][k][h], phi[i][N2-3][k][h]), phiOld); 
            error = max(error, phiOld - phi[i][N2-1][k][h]);
            
        }
        for ( i = 0; i <= (N1-1); i++)
        for ( j = 0; j <= (N2-1); j++)
        for ( k = 0; k <= (N3-1); k++)
        {
            phiOld = phi[i][j][k][0]; 
            phi[i][j][k][0] = min(max(2*phi[i][j][k][1] - phi[i][j][k][2], phi[i][j][k][2]), phiOld); 
            error = max(error, phiOld - phi[i][j][k][0]);
            
            phiOld = phi[i][j][k][N4-1];
            phi[i][j][k][N4-1] = min(max(2*phi[i][j][k][N4-2] - phi[i][j][k][N4-3], phi[i][j][k][N4-3]), phiOld); 
            error = max(error, phiOld - phi[i][j][k][N4-1]);
        }
        
//         for ( i = 0; i <= (N1-1); i++)
//         for ( j = 0; j <= (N2-1); j++)
//         {
//             phiOld = phi[i][j][0]; 
//             phi[i][j][0] = min(max(2*phi[i][j][1] - phi[i][j][2], phi[i][j][2]), phiOld); 
//             error = max(error, phiOld - phi[i][j][0]);
//             
//             phiOld = phi[i][j][N3-1];
//             phi[i][j][N3-1] = min(max(2*phi[i][j][N3-2] - phi[i][j][N3-3], phi[i][j][N3-3]), phiOld); 
//             error = max(error, phiOld - phi[i][j][N3-1]);
//         }
    }
    mexPrintf("count = %d ", count);
    return error;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //---Inside mexFunction---
    //Declarations
    double *phiValues, *xsValues, *dxValues;
    double ****phi, *****xs, *dx;
    //const int *N;
    const mwSize *N;
	int numIter;
    
    //double alpha,beta,V1,V2,error,TOL;
    //int i,j,k,l, N1,N2,N3;
    
    double a_upper,a_lower,omega_max;
    double error,TOL; 
    double *obsValues, ****obs;
    int i,j,k,h,l, N1,N2,N3,N4;
    
    //double L = 1;
    double numInfty = 1000;
    
    //Get the input
    phiValues   = mxGetPr(prhs[0]);
    xsValues    = mxGetPr(prhs[1]);
    dxValues    = mxGetPr(prhs[2]);
    
    a_upper         = (double)mxGetScalar(prhs[3]);
    a_lower         = (double)mxGetScalar(prhs[4]);
    omega_max       = (double)mxGetScalar(prhs[5]);

    
    numIter     = (int)mxGetScalar(prhs[6]);
    TOL         = (double)mxGetScalar(prhs[7]);
    obsValues       = mxGetPr(prhs[8]);
    
//     alpha       = (double)mxGetScalar(prhs[3]);
//     beta        = (double)mxGetScalar(prhs[4]);
//     V1          = (double)mxGetScalar(prhs[5]);
//     V2          = (double)mxGetScalar(prhs[6]);
    
    
    N           = mxGetDimensions(prhs[0]);
    
    N1 = N[0]; N2 = N[1]; N3 = N[2]; N4 = N[3];
    // mexPrintf("here first place");
    
    // memory allocation & value assignment
	phi   = (double ****) malloc (N1 * sizeof(double***)); 
    dx    = (double *) malloc (4 * sizeof(double)); 
    xs    = (double *****) malloc (N1 * sizeof(double****));
    obs   = (double ****) malloc ( N1 * sizeof(double***));

    // mexPrintf("here bong place");
    
	for (i=0; i<N1; i++)
    {
		phi[i]   = (double ***) malloc ( N2 * sizeof(double**));
        obs[i]   = (double ***) malloc ( N2 * sizeof(double**));
        xs[i]    = (double ****) malloc ( N2 * sizeof(double***));
        for (j=0; j<N2; j++)
        {
            phi[i][j]   = (double **) malloc ( N3 * sizeof(double*));
            obs[i][j]   = (double **) malloc ( N3 * sizeof(double*));
            xs[i][j]    = (double ***) malloc ( N3 * sizeof(double**));
            for (k=0; k<N3; k++)
            {
                phi[i][j][k] = (double *) malloc ( N4 * sizeof(double));
                obs[i][j][k] = (double *) malloc ( N4 * sizeof(double));
                xs[i][j][k]  = (double **) malloc ( N4 * sizeof(double*));
                for (h=0; h<N4; h++)
                {
                    phi[i][j][k][h] = phiValues[((h*N3+k)*N2+j)*N1+i];
                    obs[i][j][k][h] = obsValues[((h*N3+k)*N2+j)*N1+i];
                    xs[i][j][k][h]  = (double *) malloc (4 * sizeof(double));
                    for (l=0; l<4; l++)
                    {
                        xs[i][j][k][h][l] = xsValues[(((l*N4+h)*N3+k)*N2+j)*N1+i];
                    }
                }
                
            }
        }
    }
    // mexPrintf("here second place");  
    for (i=0; i<4; i++)
    {
        dx[i] = dxValues[i];
    }
    
    // run LF sweeping algorithm
    for(k=0; k<numIter; k++) 
    {
        error = sweepingLF(phi, N1,N2,N3,N4, xs, dx, a_upper,a_lower, omega_max,obs);
        // mexPrintf("here third place");
        
        mexPrintf("Error = %g at iteration %i. \n", error, k);
        mexEvalString("drawnow;");
        if (error <= TOL) {
            mexPrintf("Stopped at iteration %i. \n", k);
            break;
        }
    }
    
	// mexPrintf("here fourth place");	
  
    // send the processed phi to the output  
    for (i=0; i < N1; i++)
	for (j=0; j < N2; j++)
    for (k=0; k < N3; k++)
    for (h=0; h < N4; h++)
        //phiValues[((k*N2+j)*N1)+i] = phi[i][j][k];
          phiValues[((h*N3+k)*N2+j)*N1+i] = phi[i][j][k][h]; 
    
    // mexPrintf("here five place");
    // delete memory;
// 	for(i=0; i< N1; i++)
//     {
//         for(j=0; j<N2; j++)
//         {
//             for(k=0; k<N3; k++)
//             {
//                 free(xs[i][j][k]);
//             }
//             free(phi[i][j]);
//             free(xs[i][j]);
//         }
//         free(phi[i]);
//         free(xs[i]);
// 	}
// 	free(phi);
//     free(xs);
//     free(dx);
    
    for(i=0; i<N1; i++)
    {
        for(j=0; j<N2; j++)
        {
            for(k=0; k<N3; k++)
            {
                for(h=0; h<N4; h++)
                {
                    free(xs[i][j][k][h]);
                }
                free(phi[i][j][k]);
                free(obs[i][j][k]);
                free(xs[i][j][k]);
            }
            free(phi[i][j]);
            free(obs[i][j]);
            free(xs[i][j]);
        }
        free(phi[i]);
        free(obs[i]);
        free(xs[i]);
    }
    free(phi);
    free(obs);
    free(xs);
    free(dx);
    
}