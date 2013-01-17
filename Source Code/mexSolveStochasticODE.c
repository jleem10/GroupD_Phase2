#include "mexSolveStochasticODE.h"

// Solving the parallel system of ODEs with stochastic variables
void mexFunction(int nrhs, mxArray *prhs[], int nlhs, const mxArray *plhs[])
{

    ModelParametersStruct params;
    double *initialValues;
    size_t numInitialValues, i;
    size_t mDimNoise, nDimNoise;
    double D;
    double x[10], dx[10], xNew[10]; // To do: malloc these
    double startOfRange, endOfRange;
    double *rangeArray;
    mxArray *tmp;
    double stepSize;
    int iterations, j, k;
    double *timeOutputMatrix, *outputArray;
    double **noiseData, *noiseValues, *noise;
    
    /********************************/
    /* Input reading and validation */
    /********************************/
    
    /* Parameters */
    params.gamma_V     = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "gamma_V"));
    params.gamma_VA    = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "gamma_VA"));
    params.gamma_VH    = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "gamma_VH"));
    params.alpha_V     = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "alpha_V"));
    params.a_V1        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_V1"));
    params.a_V2        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_V2"));
    params.b_HD        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_HD"));
    params.a_R         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_R"));
    params.gamma_HV    = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "gamma_HV"));
    params.b_HF        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_HF"));
    params.b_IE        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_IE"));
    params.a_I         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_I"));
    params.b_MD        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_MD"));
    params.b_MV        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_MV"));
    params.a_M         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_M"));
    params.b_F         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_F"));
    params.c_F         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "c_F"));
    params.b_FH        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_FH"));
    params.a_F         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_F"));
    params.b_EM        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_EM"));
    params.b_EI        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_EI"));
    params.a_E         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_E"));
    params.b_PM        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_PM"));
    params.a_P         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_P"));
    params.b_A         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_A"));
    params.gamma_AV    = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "gamma_AV"));
    params.a_A         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_A"));
    params.r           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "r"));

    //Debug
    //mexPrintf("a_A: %f\n", params.a_A);
    
    /* Noise Input Values */
    noiseValues = mxGetPr(plhs[NOISE_OPTIONS]);
    nDimNoise = mxGetN(plhs[NOISE_OPTIONS]);
    mDimNoise = mxGetM(plhs[NOISE_OPTIONS]);
    
    if(mDimNoise != 2)
    {
        mexErrMsgTxt("need a 10 by 2 array with noise values\n");
        return;
    }
    
    // Allocating noiseData memory
    // Note: i is size_t
    noiseData = (double**) malloc((size_t)(mDimNoise*sizeof(double*)));
    for( i=0 ; i<mDimNoise ; i++) 
        noiseData[(int)i] = (double*)malloc((size_t)(nDimNoise*sizeof(double)));
                
    // Copying noiseData to useable format
    // Note j, k are int
    for( j=0 ; j<(int)nDimNoise ; j++)
        for( k=0 ; k<(int)mDimNoise ; k++)
            noiseData[k][j] = *(noiseValues + (size_t)k + (size_t)j*mDimNoise);
    
//      //Debug
//     mexPrintf("Noise Data\n");
//     for(i = 0; i < nDimNoise; i++)
//         mexPrintf("%f  %f  %f\n", noiseData[0][i], noiseData[3][i], noiseData[2][i]);
//     mexPrintf("\n");
    
    // Noise values to be calculated
    noise = (double*) malloc((size_t)(nDimNoise*sizeof(double)));
    
    /* Initial Values */
    initialValues = mxGetPr(plhs[INITIAL_VALUES_INDEX]);
    numInitialValues = mxGetN(plhs[INITIAL_VALUES_INDEX]);
    
    // Copying values into x
    for( i=0 ; i<numInitialValues ; i++)
    {
        x[i] = initialValues[i];
    }
    
     //Debug
//      mexPrintf("Initial Values\n");
//      for(i = 0; i < numInitialValues; i++)
//          mexPrintf("%f ", x[i]);
//      mexPrintf("\n");
    
    /* time */
    // Get the time range over which a solution should be computed
    rangeArray = mxGetPr(plhs[RANGE_INDEX]);
    if(mxGetN(plhs[RANGE_INDEX]) != 2)
    {
        mexErrMsgTxt("Time range array needs to have two element");
        return;
    }
    
    startOfRange = rangeArray[0];
    endOfRange = rangeArray[1];
    
    // Get the timestep size
    tmp = mxGetField(plhs[PROGRAM_OPTIONS], 0, "stepSize");
    if(tmp == NULL)
    {
        mexErrMsgTxt("programOptions.stepSize must be the fourth argument");
        return;
    }
    stepSize = mxGetScalar(tmp);

    // Calculating the number of Iterations 
    iterations = (int) ((endOfRange - startOfRange) / stepSize) ;
    
//      //Debug
//      mexPrintf("Iterations: %d\n", iterations);
    
    /************************/
    /* Creating Data Output */
    /************************/
    
    // Creating Time Output data
    prhs[TIME_OUTPUT_INDEX] = mxCreateDoubleMatrix((size_t)(iterations+1), 1, mxREAL);
    timeOutputMatrix = mxGetPr(prhs[TIME_OUTPUT_INDEX]);
    
    // Writing time data to output
    for( j=0 ; j<=iterations ; j++ )
    {
        timeOutputMatrix[j] = startOfRange + j*stepSize;
    }
    
    // Creating Solution output matrix
    prhs[SOLUTION_OUTPUT_START_INDEX] = mxCreateDoubleMatrix((size_t)(iterations+1), (int)numInitialValues, mxREAL);  
    outputArray = mxGetPr(prhs[SOLUTION_OUTPUT_START_INDEX]);
    
    // Saving the initial values in the output array
    for( j=0 ; j<=numInitialValues ; j++ )
        outputArray[j*(iterations+1)] = x[j];
    
    /**********************/
    /* Solving the system */
    /**********************/
    
    for( j=1 ; j<=iterations ; j++) // The main loop
    {
        
        /* Calculating changes in the quantities*/
        D = 1 - H - R - I; // D is the proportion of dead cells
        
        dV  = (params.gamma_V * I) - (params.gamma_VA * S * A * V) - (params.gamma_VH * H * V) - (params.alpha_V * V) - (params.a_V1 * V)/(1 + params.a_V2 * V);  
        dH  = (params.b_HD * D)*(H + R) + (params.a_R * R) - (params.gamma_HV * V * H) - (params.b_HF * F * H);
        dI  = (params.gamma_HV * V * H) - (params.b_IE * I * E) - (params.a_I * I);
        dM  = (params.b_MD * D + params.b_MV * V)*(1 - M) - (params.a_M * M);         
        dF  = (params.b_F * M) + (params.c_F * I) - (params.b_FH * H * F) - (params.a_F * F);
        dR = (params.b_HF * F * H) - (params.a_R * R);
        dE  = (params.b_EM * M * E) - (params.b_EI * I * E) + (params.a_E)*(1 - E);
        dP  = (params.b_PM * M * P) + (params.a_P)*(1 - P);
        dA  = (params.b_A * P) - (params.gamma_AV * S * A * V) - (params.a_A * A);
        dS  = (params .r*P)*(1 - S);
        
        /* Solving the ODE */
        
        // Noise calculation
        NoiseFunction(noiseData, (int)nDimNoise, x, noise);
        
        // ODE solution
        for( k=0 ; k<(int)numInitialValues ; k++ )
        {
            xNew[k] = dx[k]*stepSize + x[k] + noise[k]*GaussianVariable()*sqrt(stepSize);
            
            // Check for negative values in xNew, if negative --> set to 0
            // This can happen due to the noise
            if(xNew[k] < 0)
                xNew[k] = 0;
            
            /***************/ 
            /* Data output */
            /***************/
            outputArray[k*(iterations+1) + j] = xNew[k];
            
            // Preparing for next step
            x[k] = xNew[k];
            
        }        
    }   // End of main for loop
    
}

// Generates a gaussian distributed noise variable between 0 and 1
// The Box-Mueller Method is used
double GaussianVariable()
{
    double randVar;
    randVar = 0;
    
    
    //Box-Mueller equation for a Gaussian random var [0;1]
    randVar = sqrt( (-2) * log(rand0())) * cos(2 * PI * rand0());
    
    return randVar;

}

// Implementing a noise function after G(x) = ax^b
void NoiseFunction(double **noiseData, int length, double *x, double *noise)
{
    int i;
    
    for( i=0 ; i<length ; i++)
    {
        noise[i] = noiseData[0][i] * pow(x[i], noiseData[1][i]);
    }
    
    return;
}

// Random number generation in interval [0;1[
// Note: this is uniformly distributed
double rand0()
{
    double randNum;
    
    // Random seed for Gaussian Variable generation
    srand(time(NULL));
    
    randNum = (rand()%1001) / 1000.0;

    return randNum;
}