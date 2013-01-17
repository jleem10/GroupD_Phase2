#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "matrix.h"

#include "mex.h"

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

// Similar definition as above for the change in the Variables
#define dV dx[0]
#define dH dx[1]
#define dI dx[2]
#define dM dx[3]
#define dF dx[4]
#define dR dx[5]
#define dE dx[6]
#define dP dx[7]
#define dA dx[8]
#define dS dx[9]

#define PI (3.141592653589793)

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
const int NOISE_OPTIONS = 4;

// Indices of outputs in plhs[]
const int TIME_OUTPUT_INDEX = 0;
const int SOLUTION_OUTPUT_START_INDEX = 1;

// Samples a normally distribution in boundaries [0;1]
double GaussianVariable();

// Outputs the noise function G(x) = ax^b
void NoiseFunction(double **noiseData, int length, double *x, double *noise);