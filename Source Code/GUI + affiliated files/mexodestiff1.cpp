// This mex file out-sources the computation of a solution to the system of 
// differential equations being modelled to essentially be completely 
// solved in C++. 
// The INPUT is of the form
// mexodestiff(ModelParameters, [start_of_range, end_of_range], [x_0 y_0], options);
// ModelParameters is a Matlab structure containing the fields of the model 
// as defined and initialised in setInitialValue.m.

/*************************************************************************/
// start/end_of_range are the start and end values of the range to compute 
// the solution over.
// x_0 and y_0 are the initial values of the functions in the model, the 
// number of values in the vector define the number equations in the model
// and options is a Matlab structure containing options.setpSize, which 
// defines the intial step size to use for adaptive solvers.
//
// The OUTPUT is [times, solutions]
// where "times" is a vector of times
// and "solutions" is a matrix where each column is a solution for each 
// equation where the points correspond to the times in "times"
//
// % Copyright (C) 2013 Will Smith et al.
// %
// % -------------------------------------------------------------------------
// %
// %   This code is based on work by 2013 Paul Taylor et al. on
// %   on reproducing the results of Hancioglu et al.'s Paper 'A Dynamical
// %   Model of Human Immune Response to Influenza A Virus Infection'
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
// This code has used the Boost library v.1.53.0 beta 1:
// http://www.boost.org/users/history/version_1_53_0.html

#include "mex.h"

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

#include "mexodestiff1.hpp"

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
	params.b_PB        = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_PB"));
	params.a_P         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_P"));
	params.b_A         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "b_A"));
	params.gamma_AV    = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "gamma_AV"));
	params.a_A         = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "a_A"));
	params.r           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "r"));
    params.kPT           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "kPT"));
    params.kDT           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "kDT"));
    params.kPB           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "kPB"));
    params.kDB           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "kDB"));
    params.d_max           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "d_max"));
    params.t_peak           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "t_peak"));
    params.t_half           = mxGetScalar(mxGetField(plhs[MODEL_PARAMETER_INDEX], 0, "t_half"));
    
    // Get the initial values vector from Matlab
    double * initialValues = mxGetPr(plhs[INITIAL_VALUES_INDEX]);
    size_t numInitialValues = mxGetN(plhs[INITIAL_VALUES_INDEX]);
    std::vector< double > initialValuesVector(initialValues, initialValues + numInitialValues);

//     //DEBUG
//     mexPrintf("Initial Values\n");
//     for(size_t i = 0; i < numInitialValues; i++)
//         mexPrintf("%f ", initialValues[i]);
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
    
    // Get the initial step size (or step size of a non-adaptive solver is being used)
    mxArray * tmp = mxGetField(plhs[PROGRAM_OPTIONS], 0, "stepSize");
    if(tmp == NULL)
    {
        mexErrMsgTxt("programOptions.stepSize must be the fourth argument");
        return;
    }
    double initialStepSize = mxGetScalar(tmp);

    // Copy the intial values into a vector for odeint
    VectorType x(initialValuesVector.size());
    for(size_t i = 0; i < initialValuesVector.size(); i++)
    {
        x[i] = initialValuesVector[i];
    }
    
    // These hold the solutions and times
    std::vector< VectorType > xVec;
    std::vector< double >  times;
    
    // Use the rosenbrock4 solver (recommended for stiff systems) with the systems defined in StiffSystem and Jacobian defined in StiffSystemJacobian
    // The arguments to make_dense_output set the absolute and relative error used in the computation by the solver
    // x contains the initial values
    // PushBackStateAndTime is used as a way to store the solution and time as the solver progresses
    size_t steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-3 ), make_pair(StiffSystem(params), StiffSystemJacobian(params)), 
                x, startOfRange, endOfRange, initialStepSize, PushBackStateAndTime(xVec, times));
    
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
