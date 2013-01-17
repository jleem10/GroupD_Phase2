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

#include <vector>
#include <utility>
using namespace std;

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

// Similar to the above definitions, but just the index for each function
// Helps to make the Jacobian more readable
#define Vi 0
#define Hi 1
#define Ii 2
#define Mi 3
#define Fi 4
#define Ri 5
#define Ei 6
#define Pi 7
#define Ai 8
#define Si 9

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
} ModelParametersStruct;

// Indices of input parameters in prhs[]
const int MODEL_PARAMETER_INDEX = 0;
const int RANGE_INDEX = 1;
const int INITIAL_VALUES_INDEX = 2;
const int PROGRAM_OPTIONS = 3;

// Indices of outputs in plhs[]
const int TIME_OUTPUT_INDEX = 0;
const int SOLUTION_OUTPUT_START_INDEX = 1;

// This definition is required to use odeint's stiff solvers
typedef boost::numeric::ublas::vector< double > VectorType;
typedef boost::numeric::ublas::matrix< double > MatrixType;

// odeint requires this structure to define the system of differential equations being modelled
// Using a struct also allows the parameters to be stored as a local variable rather than being global
struct StiffSystem
{
    ModelParametersStruct params;
    
    StiffSystem(ModelParametersStruct params) : params(params) {}
    
    // For odeint, the () operator is overloaded with the system of differential equations
    // Because the preprocessor has been used to define the functions (V, H, I etc.) in terms of a vector "x" the input argument must have the name x
    void operator()( const VectorType &x, VectorType &dxdy, double /* t */ )
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
        dxdy[9]  = (params .r*P)*(1 - S);
    }
};

// odeint requires the Jacobian of the stiff system of differential equations to be defined as a structure in the following manner
struct StiffSystemJacobian
{
    ModelParametersStruct params;
    StiffSystemJacobian(ModelParametersStruct params) : params(params) {}
    
    // odeint requires that the () operator be overloaded to write the Jacobian matrix out
    void operator()(const VectorType & x , MatrixType &J , const double & /* t */ , VectorType &dfdt)
    {
        double D = 1 - H - R - I;
        
        J(Vi, Vi) = -params.gamma_VA * S * A - params.gamma_VH * H - params.alpha_V - (params.a_V1/(1 + params.a_V2 * V));
        J(Vi, Hi) = -params.gamma_VH * V;
        J(Vi, Ii) = -params.gamma_V;
        J(Vi, Mi) = 0;
        J(Vi, Fi) = 0;
        J(Vi, Ri) = 0;
        J(Vi, Ei) = 0;
        J(Vi, Pi) = 0;
        J(Vi, Ai) = -params.gamma_VA * S * V;
        J(Vi, Si) = -params.gamma_VA * A * V;
        
        J(Hi, Vi) = -params.gamma_HV * V;
        J(Hi, Hi) = -params.b_HD * D - params.gamma_HV * V;
        J(Hi, Ii) = 0;
        J(Hi, Mi) = 0;
        J(Hi, Fi) = -params.b_HF * H;
        J(Hi, Ri) = params.b_HD * D + params.a_R;
        J(Hi, Ei) = 0;
        J(Hi, Pi) = 0;
        J(Hi, Ai) = 0;
        J(Hi, Si) = 0;

        J(Ii, Vi) = params.gamma_HV * H;
        J(Ii, Hi) = params.gamma_HV * V;
        J(Ii, Ii) = -params.b_IE * E - params.a_I;
        J(Ii, Mi) = 0;
        J(Ii, Fi) = 0;
        J(Ii, Ri) = 0;
        J(Ii, Ei) = -params.b_IE * I;
        J(Ii, Pi) = 0;
        J(Ii, Ai) = 0;
        J(Ii, Si) = 0;
        
        J(Mi, Vi) = params.b_MV - params.b_MV * M;
        J(Mi, Hi) = 0;
        J(Mi, Ii) = 0;
        J(Mi, Mi) = params.b_MD * D - params.b_MV * V - params.a_M;
        J(Mi, Fi) = 0;
        J(Mi, Ri) = 0;
        J(Mi, Ei) = 0;
        J(Mi, Pi) = 0;
        J(Mi, Ai) = 0;
        J(Mi, Si) = 0;
        
        J(Fi, Vi) = 0;
        J(Fi, Hi) = -params.b_FH * F;
        J(Fi, Ii) = params.c_F;
        J(Fi, Mi) = params.b_F;
        J(Fi, Fi) = -params.b_FH * H - params.a_F;
        J(Fi, Ri) = 0;
        J(Fi, Ei) = 0;
        J(Fi, Pi) = 0;
        J(Fi, Ai) = 0;
        J(Fi, Si) = 0;
        
        J(Ri, Vi) = 0;
        J(Ri, Hi) = params.b_HF * F;
        J(Ri, Ii) = 0;
        J(Ri, Mi) = 0;
        J(Ri, Fi) = params.b_HF * H;
        J(Ri, Ri) = -params.a_R;
        J(Ri, Ei) = 0;
        J(Ri, Pi) = 0;
        J(Ri, Ai) = 0;
        J(Ri, Si) = 0;
        
        J(Ei, Vi) = 0;
        J(Ei, Hi) = 0;
        J(Ei, Ii) = -params.b_EI * E;
        J(Ei, Mi) = params.b_EM * E;
        J(Ei, Fi) = 0;
        J(Ei, Ri) = 0;
        J(Ei, Ei) = params.b_EM * M - params.b_EI * I - params.a_E;
        J(Ei, Pi) = 0;
        J(Ei, Ai) = 0;
        J(Ei, Si) = 0;
        
        J(Pi, Vi) = 0;
        J(Pi, Hi) = 0;
        J(Pi, Ii) = 0;
        J(Pi, Mi) = params.b_PM * P;
        J(Pi, Fi) = 0;
        J(Pi, Ri) = 0;
        J(Pi, Ei) = 0;
        J(Pi, Pi) = params.b_PM * M + params.a_P;
        J(Pi, Ai) = 0;
        J(Pi, Si) = 0;
        
        J(Ai, Vi) = -params.gamma_AV * S * A;
        J(Ai, Hi) = 0;
        J(Ai, Ii) = 0;
        J(Ai, Mi) = 0;
        J(Ai, Fi) = 0;
        J(Ai, Ri) = 0;
        J(Ai, Ei) = 0;
        J(Ai, Pi) = params.b_A;
        J(Ai, Ai) = -params.gamma_AV * S * V - params.a_A;
        J(Ai, Si) = -params.gamma_AV * A * V;
        
        J(Si, Vi) = 0;
        J(Si, Hi) = 0;
        J(Si, Ii) = 0;
        J(Si, Mi) = 0;
        J(Si, Fi) = 0;
        J(Si, Ri) = 0;
        J(Si, Ei) = 0;
        J(Si, Pi) = params.r - S;
        J(Si, Ai) = 0;
        J(Si, Si) = -params.r * P;
        
        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
        dfdt[2] = 0.0;
        dfdt[3] = 0.0;
        dfdt[4] = 0.0;
        dfdt[5] = 0.0;
        dfdt[6] = 0.0;
        dfdt[7] = 0.0;
        dfdt[8] = 0.0;
        dfdt[9] = 0.0;
    }
};

// This struture is used by odeint to keep track of the solution at each time step, so it can be written out to Matlab
struct PushBackStateAndTime
{
    // "states" contains the solutions at the time in the corresponding index of "times"
    std::vector< VectorType >& states;
    std::vector< double >& times;
    
    PushBackStateAndTime(std::vector< VectorType >& states, std::vector< double >& times)
        : states(states), times(times) {}
    
    // odeint uses this overloaded operator to pass the computed solution at each time step
    void operator() (const VectorType& x, double t)
    {
        states.push_back(x);
        times.push_back(t);
    }
};
