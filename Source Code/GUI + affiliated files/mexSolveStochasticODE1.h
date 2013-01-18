// % -------------------------------------------------------------------------
// %
// %   A Program to solve a stochastic implementation of the system of ODEs
// %   first proposed by Hancioglu et al. in the paper 'A Dynamical Model
// %   of Human Immune Response to Influenza A Virus Infection'
// %   Copyright (C) 2013 Will Smith et al.
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
// %
// %   This code is based on work by 2013 Paul Taylor et al. on
// %   on reproducing the results of Hancioglu et al.'s Paper 'A Dynamical
// %   Model of Human Immune Response to Influenza A Virus Infection'


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
// % TV         T cell clone
// % Bi         Activated B cell clone
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
#define TV x[10]
#define Bi x[11]

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
#define dTV dx[10]
#define dBi dx[11]

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
	double b_PB;
	double a_P;
	double b_A;
	double gamma_AV;
	double a_A;
	double r;
    double kPT;
    double kDT;
    double kPB;
    double kDB;
    double d_max;
    double t_peak;
    double t_half;
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

// Samples a normal distribution N(0;1)
double GaussianVariable();

// Outputs the noise function G(x) = ax^b
void NoiseFunction(double **noiseData, int length, double *x, double *noise);

// Outputs a uniformly distributed random number between 0 and 1
// Note: there is a bias if RAND_MAX >> 1000 is false
double rand0();