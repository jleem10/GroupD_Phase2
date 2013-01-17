%SETINITIALVALUE sets the value of the variables in the initialisation tab
%to a new value specified by the user.
%
% handles.initial_parameter is the number of the selected variable name on
% the popupmenu. This function associates that number with the correct
% element of the vector, InitialParameters, or the struct, ModelParameters.
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

function handles = setInitialValue(handles)

if (handles.initial_parameter>=1 && handles.initial_parameter<=10)
    handles.InitialConditions(handles.initial_parameter) = handles.new_value;
else
   switch handles.initial_parameter
       case 11
            handles.ModelParameters.gamma_V = handles.new_value;
       case 12
            handles.ModelParameters.gamma_VA = handles.new_value;
       case 13
            handles.ModelParameters.gamma_VH = handles.new_value;
       case 14
            handles.ModelParameters.alpha_V = handles.new_value;
       case 15
            handles.ModelParameters.a_V1 = handles.new_value;
       case 16
            handles.ModelParameters.a_V2 = handles.new_value;
       case 17
            handles.ModelParameters.b_HD = handles.new_value;
       case 18
            handles.ModelParameters.a_R = handles.new_value;
       case 19
            handles.ModelParameters.gamma_HV = handles.new_value;
       case 20
            handles.ModelParameters.b_HF = handles.new_value;
       case 21
            handles.ModelParameters.b_IE = handles.new_value;
       case 22
            handles.ModelParameters.a_I = handles.new_value;
       case 23
            handles.ModelParameters.b_MD = handles.new_value;
       case 24
            handles.ModelParameters.b_MV = handles.new_value;
       case 25
            handles.ModelParameters.a_M = handles.new_value;
       case 26
            handles.ModelParameters.b_F = handles.new_value;
       case 27
            handles.ModelParameters.c_F = handles.new_value;
       case 28
            handles.ModelParameters.b_FH = handles.new_value;
       case 29
            handles.ModelParameters.a_F = handles.new_value;
       case 30
            handles.ModelParameters.b_EM = handles.new_value;
       case 31
            handles.ModelParameters.b_EI = handles.new_value;
       case 32
            handles.ModelParameters.a_E = handles.new_value;
       case 33
            handles.ModelParameters.b_PM = handles.new_value;
       case 34
            handles.ModelParameters.a_P = handles.new_value;
       case 35
            handles.ModelParameters.b_A = handles.new_value;
       case 36
            handles.ModelParameters.gamma_AV = handles.new_value;
       case 37
            handles.ModelParameters.a_A = handles.new_value;
       case 38
         handles.ModelParameters.r = handles.new_value;
       case 39
         handles.ModelParameters.d_max = handles.new_value;
       case 40
         handles.ModelParameters.t_peak = handles.new_value;
       case 41
         handles.ModelParameters.t_half = handles.new_value;
   end
end