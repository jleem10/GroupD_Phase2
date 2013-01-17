%
%   -----------------------------------------------------------------------
%
%   A Program to implement and reproduce the results of Hancioglu et al.'s Paper
%   'A Dynamical Model of Human Immune Response to Influenza A Virus Infection'
%   Copyright (C) 2013  Paul Taylor et al.
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [ Solution ] = SolveODE_simple( InitialConditions )
%Function to Solve the Differential Equations Governing Model Variable
%Populations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Differential Variables:
%
% Y(1) = V          Viral Load Per Epithelial Cell
% Y(2) = H          Proportion of Healthy Cells
% Y(3) = I          Proportion of Infected Cells
% Y(4) = M          Activated Antigen-Presenting Cells Per Homeostatic Level
% Y(5) = F          Interferons Per Homeostatic Level Of Macrophages
% Y(6) = R          Proportion of Resistant Cells
% Y(7) = E          Effector Cells per Homeostatic Level 
% Y(8) = P          Plasma Cells Per Homeostatic Level 
% Y(9) = A          Antibodies Per Homeostatic Level   
% Y(10)= S          Antigenic Distance
% 
% D                 Proportion of Dead Cells
%
% D is not listed as a Variable in Y as it is solely determined 
% % by (1-H-R-I) and so can be directly inferred during afterwards
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ModelParameters is a Struct with elements
% 
% ModelParameters.gamma_V
% ModelParameters.gamma_VA
% ModelParameters.gamma_VH
% ModelParameters.alpha_V
% ModelParameters.a_V1
% ModelParameters.a_V2
% ModelParameters.b_HD
% ModelParameters.a_R
% ModelParameters.gamma_HV
% ModelParameters.b_HF
% ModelParameters.b_IE
% ModelParameters.a_I
% ModelParameters.b_MD
% ModelParameters.b_MV
% ModelParameters.a_M
% ModelParameters.b_F
% ModelParameters.c_F
% ModelParameters.b_FH
% ModelParameters.a_F
% ModelParameters.b_EM
% ModelParameters.b_EI
% ModelParameters.a_E
% ModelParameters.b_PM
% ModelParameters.a_P
% ModelParameters.b_A
% ModelParameters.gamma_AV
% ModelParameters.a_A
% ModelParameters.r
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ProgramParameters is a Struct with elements
% 
% ProgramParameters.t_start
% ProgramParameters.t_end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ProgramParameters;

[Solution.t,Solution.Y] = ode15s(@DifferentialEquations,[ProgramParameters.t_start,ProgramParameters.t_end],InitialConditions);

Solution.Y(:,11) = 1 - Solution.Y(:,2) - Solution.Y(:,6) - Solution.Y(:,3);


plot(Solution.Y(:,11),Solution.Y(:,1));

end

function [ dYdt ] = DifferentialEquations( t , Y )
%Differential Equations Governing the Model Populations

global ModelParameters;

V = Y(1);
H = Y(2);
I = Y(3);
M = Y(4);
F = Y(5);
R = Y(6);
E = Y(7);
P = Y(8);
A = Y(9);
S = Y(10);
D = 1 - H - R - I;

dYdt(1,1)  = (ModelParameters.gamma_V * I) - (ModelParameters.gamma_VA * S * A * V) - (ModelParameters.gamma_VH * H * V) - (ModelParameters.alpha_V * V) - (ModelParameters.a_V1 * V)/(1 + ModelParameters.a_V2 * V);
dYdt(2,1)  = (ModelParameters.b_HD * D)*(H + R) + (ModelParameters.a_R * R) - (ModelParameters.gamma_HV * V * H) - (ModelParameters.b_HF * F * H);
dYdt(3,1)  = (ModelParameters.gamma_HV * V * H) - (ModelParameters.b_IE * I * E) - (ModelParameters.a_I * I);
dYdt(4,1)  = (ModelParameters.b_MD * D + ModelParameters.b_MV * V)*(1 - M) - (ModelParameters.a_M * M);
dYdt(5,1)  = (ModelParameters.b_F * M) + (ModelParameters.c_F * I) - (ModelParameters.b_FH * H * F) - (ModelParameters.a_F * F);
dYdt(6,1)  = (ModelParameters.b_HF * F * H) - (ModelParameters.a_R * R);
dYdt(7,1)  = (ModelParameters.b_EM * M * E) - (ModelParameters.b_EI * I * E) + (ModelParameters.a_E)*(1 - E);
dYdt(8,1)  = (ModelParameters.b_PM * M * P) + (ModelParameters.a_P)*(1 - P);
dYdt(9,1)  = (ModelParameters.b_A * P) - (ModelParameters.gamma_AV * S * A * V) - (ModelParameters.a_A * A);
dYdt(10,1) = (ModelParameters.r*P)*(1 - S);

end








