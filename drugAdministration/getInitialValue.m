%GETINITIALVALUE returns the current value of a parameter held within the
% ModelParameters struct or the InitialConditions vector. These values can
% be set by the user in the GUI.
%
% -------------------------------------------------------------------------
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

function current_value = getInitialValue(handles,parameter)

if (parameter>=1 && parameter<=10)
    current_value = handles.InitialConditions(parameter);
else
   switch parameter
       case 11
            current_value = handles.ModelParameters.gamma_V;
       case 12
            current_value = handles.ModelParameters.gamma_VA;
       case 13
            current_value = handles.ModelParameters.gamma_VH;
       case 14
            current_value = handles.ModelParameters.alpha_V;
       case 15
            current_value = handles.ModelParameters.a_V1;
       case 16
            current_value = handles.ModelParameters.a_V2;
       case 17
            current_value = handles.ModelParameters.b_HD;
       case 18
            current_value = handles.ModelParameters.a_R;
       case 19
            current_value = handles.ModelParameters.gamma_HV;
       case 20
            current_value = handles.ModelParameters.b_HF;
       case 21
            current_value = handles.ModelParameters.b_IE;
       case 22
            current_value = handles.ModelParameters.a_I;
       case 23
            current_value = handles.ModelParameters.b_MD;
       case 24
            current_value = handles.ModelParameters.b_MV;
       case 25
            current_value = handles.ModelParameters.a_M;
       case 26
            current_value = handles.ModelParameters.b_F;
       case 27
            current_value = handles.ModelParameters.c_F;
       case 28
            current_value = handles.ModelParameters.b_FH;
       case 29
            current_value = handles.ModelParameters.a_F;
       case 30
            current_value = handles.ModelParameters.b_EM;
       case 31
            current_value = handles.ModelParameters.b_EI;
       case 32
            current_value = handles.ModelParameters.a_E;
       case 33
            current_value = handles.ModelParameters.b_PM;
       case 34
            current_value = handles.ModelParameters.a_P;
       case 35
            current_value = handles.ModelParameters.b_A;
       case 36
            current_value = handles.ModelParameters.gamma_AV;
       case 37
            current_value = handles.ModelParameters.a_A;
       case 38
            current_value = handles.ModelParameters.r;
       case 39
            current_value = handles.ModelParameters.d_max;
       case 40
            current_value = handles.ModelParameters.t_peak;
       case 41
            current_value = handles.ModelParameters.t_half;
       case 42
            current_value = handles.ModelParameters.t_on;
   end
end