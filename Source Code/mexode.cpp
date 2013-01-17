// % -------------------------------------------------------------------------
// %
// %   A Program to implement and reproduce the results of Hancioglu et al.'s Paper
// %   'A Dynamical Model of Human Immune Response to Influenza A Virus Infection'
// %   Copyright (C) 2013  Paul Taylor et al.
// % 
// %   This program is free software: you can redistribute it and/or modify
// %   it under the terms of the GNU General Public License as published by
// %   the Free Software Foundation, either version 3 of the License, or
// %   (at your option) any later version.
// % 
// %   This program is distributed in the hope that it will be useful,
// %   but WITHOUT ANY WARRANTY; without even the implied warranty of
// %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// %   GNU General Public License for more details.
// % 
// %   You should have received a copy of the GNU General Public License
// %   along with this program.  If not, see <http://www.gnu.org/licenses/>.


// This code uses the ODEINT library http://www.odeint.com 
// which is published under the boost license http://www.boost.org/LICENSE_1_0.txt

// This mex file out-sources the computation of a solution to the system of differential equations being modelled to essentially be completely solved in C++
// The solvers used in this mex file are for NON-STIFF problems
// The INPUT is of the form
// mexode(ModelParameters, [start_of_range, end_of_range], [x_0 y_0], options);
// ModelParameters is a Matlab structure containing the fields of the model as defined and initialised in InitialiseModelParameters.m
// start/end_of_range are the start and end values of the range to compute the solution over
// x_0 and y_0 are the initial values of the functions in the model, the number of values in the vector define the number equations in the model
// and options is a Matlab structure containing options.setpSize, which defines the intial step size to use for adaptive solvers
// The OUTPUT is [times, solutions]
// where "times" is a vector of times
// and "solutions" is a matrix where each column is a solution for each equation where the points correspond to the times in "times"

#include "mex.h"

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

#include <vector>
typedef std::vector< double > StateType;

// Define elements of the solution matrix using the notation
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %
// % Differential Variables:
// %
// % V          Viral Load Per Epithelial Cell
// % H          Proportion of Healthy Cells
// % I          Proportion of Infected Cells
// % M          Activated Antigen-Presenting Cells Per Homeostatic Level
// % F          Interferons Per Homeostatic Level Of Macrophages
// % R          Proportion of Resistant Cells
// % E          Effector Cells per Homeostatic Level 
// % P          Plasma Cells Per Homeostatic Level 
// % A          Antibodies Per Homeostatic Level   
// % S          Antigenic Distance
// % 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define V x[0]
#define H x[1]
#define I x[2]
#define M x[3]
#define F x[4]
#define R x[5]
#define E x[6]
#define P x[7]
#define A x[8]
#define S x[9]
#define X x[10]          // %----CHANGE----%

// Indices of input parameters in prhs[]
const int MODEL_PARAMETER_INDEX = 0;
const int RANGE_INDEX = 1;
const int INITIAL_VALUES_INDEX = 2;
const int PROGRAM_OPTIONS = 3;

// Indices of outputs in plhs[]
const int TIME_OUTPUT_INDEX = 0;
const int SOLUTION_OUTPUT_START_INDEX = 1;

// Struct to organise all the model parameters
typedef struct ModelParametersStruct
{
    double gamma_V;
    double gamma_VA;
	double gamma_VH;
	double alpha_V;
	double a_V1;
	double a_V2;
	double b_HD;
	double a_R;
	double gamma_HV;
	double b_HF;
	double b_IE;
	double a_I;
	double b_MD;
	double b_MV;
	double a_M;
	double b_F;
	double c_F;
	double b_FH;
	double a_F;
	double b_EM;
	double b_EI;
	double a_E;
	double b_PM;
	double a_P;
	double b_A;
	double gamma_AV;
	double a_A;
	double r;
    double c_1;          // %----CHANGE----%
    double c_2;          // %----CHANGE----%
    double c_3;          // %----CHANGE----%
} ModelParametersStruct;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %
// % Differential Variables:
// %
// % V          Viral Load Per Epithelial Cell
// % H          Proportion of Healthy Cells
// % I          Proportion of Infected Cells
// % M          Activated Antigen-Presenting Cells Per Homeostatic Level
// % F          Interferons Per Homeostatic Level Of Macrophages
// % R          Proportion of Resistant Cells
// % E          Effector Cells per Homeostatic Level 
// % P          Plasma Cells Per Homeostatic Level 
// % A          Antibodies Per Homeostatic Level   
// % S          Antigenic Distance
// % 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// odeint requires this class to define the system of differential equations being modelled
// Using a class also allows the parameters to be stored as a local variable rather than being global
class ODEs
{
    // Store the parameters for the model
    ModelParametersStruct params;
    
    public:
        ODEs(ModelParametersStruct modelParameters) : params(modelParameters) {} 
    
    // For odeint, the () operator is overloaded with the system of differential equations
    // Because the preprocessor has been used to define the functions (V, H, I etc.) in terms of a vector "x" the input argument must have the name x
    void operator() ( const StateType &x , StateType &dxdy , const double /* t */ )
    {
        double D = 1 - H - R - I; // D is the proportion of dead cells
        
        dxdy[0]  = (params.gamma_V * I) - (params.gamma_VA * S * A * V) - (params.gamma_VH * H * V) - (params.alpha_V * V) - (params.a_V1 * V)/(1 + params.a_V2 * V);  
        dxdy[1]  = (params.b_HD * D)*(H + R) + (params.a_R * R) - (params.gamma_HV * V * H) - (params.b_HF * F * H);
        dxdy[2]  = (params.gamma_HV * V * H) - (params.b_IE * I * E) - (params.a_I * I);
        dxdy[3]  = (params.b_MD * D + params.b_MV * V)*(1 - M) - (params.a_M * M);         
        dxdy[4]  = (params.b_F * M) + (params.c_F * I) - (params.b_FH * H * F) - (params.a_F * F);
        dxdy[5]  = (params.b_HF * F * H) - (params.a_R * R);
        dxdy[6]  = (params.b_EM * M * E) - (params.b_EI * I * E) + (params.a_E)*(1 - E);
        dxdy[7]  = (params.b_PM * M * P) + (params.a_P)*(1 - P);
        dxdy[8]  = (params.b_A * P) - (params.gamma_AV * S * A * V) - (params.a_A * A);
        dxdy[9]  = (params.r*P)*(1 - S);
        dxdy[10] = (params.c_1 * I * X) - (params.c_2 * X) + (params.c_3);          // %----CHANGE----%
    }
};

// This struture is used by odeint to keep track of the solution at each time step, so it can be written out to Matlab
struct PushBackStateAndTime
{
    // "states" contains the solutions at the time in the corresponding index of "times"
    std::vector< StateType >& states;
    std::vector< double >& times;
    
    PushBackStateAndTime(std::vector< StateType >& states, std::vector< double >& times)
        : states(states), times(times) {}
    
    // odeint uses this overloaded operator to pass the computed solution at each time step
    void operator() (const StateType& x, double t)
    {
        states.push_back(x);
        times.push_back(t);
    }
};

void mexFunction(int nrhs, mxArray *prhs[], int nlhs, const mxArray *plhs[])
{
    // Get input data from Matlab
    // TODO input validation
    // Store the mode parameters from their struct
    ModelParametersStruct params;
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
    params.c_1         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "c_1"));          // %----CHANGE----%
    params.c_2         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "c_2"));          // %----CHANGE----%
    params.c_3         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "c_3"));          // %----CHANGE----%
    
    // Get the initial values vector from Matlab
    double * initialValues = mxGetPr(plhs[INITIAL_VALUES_INDEX]);
    size_t numInitialValues = mxGetN(plhs[INITIAL_VALUES_INDEX]);
    std::vector< double > initialValuesVector(initialValues, initialValues + numInitialValues);
    
    // Debug
//     mexPrintf("Initial Values\n");
//     for(size_t i = 0; i < numInitialValues; i++)
//     {
//         mexPrintf("%f ", initialValues[i]);
//     }
//     mexPrintf("\n");

    // Get the range over which a solution should be computed
    double * rangeArray = mxGetPr(plhs[RANGE_INDEX]);
    if(mxGetN(plhs[RANGE_INDEX]) != 2)
    {
        mexErrMsgTxt("Range array needs to have two element");
        return;
    }
    
    double startOfRange = rangeArray[0];
    double endOfRange = rangeArray[1];
    
    mxArray * tmp = mxGetField(plhs[PROGRAM_OPTIONS], 0, "stepSize"); //This is just the initial step size for adaptive solvers
    if(tmp == NULL)
    {
        mexErrMsgTxt("programOptions.stepSize must be the fourth argument");
        return;
    }
    double initialStepSize = mxGetScalar(tmp);
    
    // This contains the initial values for odeint
    StateType x(initialValuesVector);
    
    // These hold the solutions and times
    std::vector< StateType > xVec;
    std::vector< double >  times;
    
    // Ititialise the model
    ODEs odes(params);
    
    // Different solvers...
    // runge_kutta
    // runge_kutta_dopri5
    // runge_kutta_cash_karp54
    
    // Using appropriate methods
    // integrage_const
    // integrate_adaptive
    // integrate_n_steps
    
//     runge_kutta4< StateType > stepper;
//     size_t steps = integrate_const(stepper, odes, x, startOfRange, endOfRange, initialStepSize, PushBackStateAndTime(xVec, times));

//     runge_kutta4< StateType > stepper;
//     size_t steps = integrate_adaptive(stepper, odes, x , startOfRange, endOfRange, initialStepSize, PushBackStateAndTime(xVec, times));
    
//     typedef runge_kutta_dopri5< double > ErrorStepperType;
//     size_t steps = integrate_adaptive(make_controlled< ErrorStepperType >( 1.0e-12 , 1.0e-6 ), odes, x , startOfRange, endOfRange, initialStepSize, PushBackStateAndTime(xVec, times));   
//     size_t steps = integrate_const( make_dense_output< runge_kutta_dopri5< StateType > >( 1.0e-6 , 1.0e-3 ), odes, x , startOfRange, endOfRange, initialStepSize, PushBackStateAndTime(xVec, times));
//     size_t steps = integrate(odes, x, startOfRange, endOfRange, initialStepSize, PushBackStateAndTime(xVec, times));

//     typedef runge_kutta_cash_karp54< StateType > ErrorStepperType;
//     typedef controlled_runge_kutta< ErrorStepperType > ControlledStepperType;
//     ControlledStepperType controlledStepper;
//     size_t steps = integrate_n_steps( controlledStepper, odes, x, startOfRange, endOfRange, 1000, PushBackStateAndTime(xVec, times));
    
    // Perform the integration
    // the numerical arguments to make_controlled specify the absolute and relative error respectively
    // odes is the system of equations
    // x contains the initial values
    // PushBackStateAndTime is used to write the solution and time step as the solution progresses
    typedef runge_kutta_cash_karp54< StateType > ErrorStepperType;
    size_t steps = integrate_adaptive(make_controlled< ErrorStepperType >( 1.0e-6, 1.0e-3 ), odes, x , startOfRange, endOfRange, initialStepSize, PushBackStateAndTime(xVec, times));   
    
    // Output Data to Matlab
    // TODO check for memory errors
    // Output the times in a vector
    prhs[TIME_OUTPUT_INDEX] = mxCreateDoubleMatrix(times.size(), 1, mxREAL);
    
    double *timeOutputMatrix = mxGetPr(prhs[TIME_OUTPUT_INDEX]);
    for(size_t i = 0; i < times.size(); i++)
    {
        timeOutputMatrix[i] = times[i];
    }
    
    // Output the solutions in a matrix in which each column is a solution to a single equation in the system
    prhs[SOLUTION_OUTPUT_START_INDEX] = mxCreateDoubleMatrix(xVec.size(), numInitialValues, mxREAL);  
    double *outputArray = mxGetPr(prhs[SOLUTION_OUTPUT_START_INDEX]);
    
    for(int i = 0; i < numInitialValues; i++)
    {
        for(size_t j = 0; j < xVec.size(); j++)
        {
            outputArray[i*xVec.size() + j] = xVec[j][i];
        }
    }
}
