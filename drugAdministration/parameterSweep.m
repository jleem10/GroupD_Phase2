%PARAMETERSWEEP performs nsweeps iterations of the ODE solver to loop over
% all the required values of the sweep variable.
%
% A parfor loop has been used to improve the efficiency of the loop, and a
% MEX function is available as an alternative solver to MATLAB's ode15s.
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

function [Output LegendString SweepVariableName, SweepValues] = parameterSweep(InitialConditions, ProgramParameters, ModelParameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get The Name of the Parameter Being Swept Over %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ProgramParameters.sweepvar
   case 1
        SweepVariableName = 'InitialConditions(1)'; % Variable Name (Used In-Function)
        Name = 'V(0)';                              % For Use In Legends
   case 2
        SweepVariableName = 'InitialConditions(2)';
        Name = 'H(0)'; 
   case 3
        SweepVariableName = 'InitialConditions(3)';
        Name = 'I(0)';
   case 4
        SweepVariableName = 'InitialConditions(4)';
        Name = 'M(0)';
   case 5
        SweepVariableName = 'InitialConditions(5)';
        Name = 'F(0)';
   case 6
        SweepVariableName = 'InitialConditions(6)';
        Name = 'R(0)';
   case 7
        SweepVariableName = 'InitialConditions(7)';
        Name = 'E(0)';
   case 8
        SweepVariableName = 'InitialConditions(8)';
        Name = 'P(0)';
   case 9
        SweepVariableName = 'InitialConditions(9)';
        Name = 'A(0)';
   case 10
        SweepVariableName = 'InitialConditions(10)';
        Name = 'S(0)';
   case 11
        SweepVariableName = 'gamma_V';
   case 12
        SweepVariableName = 'gamma_VA';
   case 13
        SweepVariableName = 'gamma_VH';
   case 14
        SweepVariableName = 'alpha_V';
   case 15
        SweepVariableName = 'a_V1';
   case 16
        SweepVariableName = 'a_V2';
   case 17
        SweepVariableName = 'b_HD';
   case 18
        SweepVariableName = 'a_R';
   case 19
        SweepVariableName = 'gamma_HV';
   case 20
        SweepVariableName = 'b_HF';
   case 21
        SweepVariableName = 'b_IE';
   case 22
        SweepVariableName = 'a_I';
   case 23
        SweepVariableName = 'b_MD';
   case 24
        SweepVariableName = 'b_MV';
   case 25
        SweepVariableName = 'a_M';
   case 26
        SweepVariableName = 'b_F';
   case 27
        SweepVariableName = 'c_F';
   case 28
        SweepVariableName = 'b_FH';
   case 29
        SweepVariableName = 'a_F';
   case 30
        SweepVariableName = 'b_EM';
   case 31
        SweepVariableName = 'b_EI';
   case 32
        SweepVariableName = 'a_E';
   case 33
        SweepVariableName = 'b_PM';
   case 34
        SweepVariableName = 'a_P';
   case 35
        SweepVariableName = 'b_A';
   case 36
        SweepVariableName = 'gamma_AV';
   case 37
        SweepVariableName = 'a_A';
   case 38
        SweepVariableName = 'r';
   case 39
        SweepVariableName = 'd_max';
   case 40
        SweepVariableName = 't_peak';
   case 41
        SweepVariableName = 't_half';
   case 42
        SweepVariableName = 't_on';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up Cells and Arrays for Simulations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ProgramParameters.logsweep == 0
    SweepValues = linspace(ProgramParameters.min,ProgramParameters.max,...
         ProgramParameters.nsweeps);
else
     SweepValues = logspace(log10(ProgramParameters.min),log10(ProgramParameters.max),...
        ProgramParameters.nsweeps);
end
    
%%% Make Cell Array To Hold Output %%% 
Output = cell(ProgramParameters.nsweeps,1);

%%% Make Cell Array for Varying Forms of ModelParameters %%%
LocalModelParameters = cell(ProgramParameters.nsweeps,1); 

%%% Produce nsweeps Copies of InitialConditions %%%
LocalInitialConditions = InitialConditions(ones(1,ProgramParameters.nsweeps),:);

%%% Make Lists of Legend Entries %%%
LegendString = cell(ProgramParameters.nsweeps,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Change Parameter Sets For Different Values of SweepVariable %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% For each set of LocalModelParameters, change the SweepVariable to one of its possible values %%%
if (ProgramParameters.sweepvar > 10)
    for Index=1:ProgramParameters.nsweeps

		%%% For The SweepVariable == A Model Parameter %%%

		%%% Make a Local Copy of ModelParameters %%%
        LocalModelParameters{Index} = ModelParameters;

		%%% Create Command in form of a String %%%
		%%% This will change the ModelParameter which corresponds to the SweepVariable to one of the SweepValues %%%
        str = sprintf('LocalModelParameters{%i}.%s = SweepValues(%i);',Index,SweepVariableName,Index);
        
		%%% Evaluate as Command %%%
		eval(str);

		%%% Create String To Use In Legend %%%
        LegendString{Index} = sprintf('%s = %.3g',SweepVariableName,SweepValues(Index));

    end
else
    for Index=1:ProgramParameters.nsweeps

		%%% For The SweepVariable == an Initial Condition %%%
		
		%%% Make a Copy of LocalModelParameters %%%
		LocalModelParameters{Index} = ModelParameters;

		%%% Change the LocalInitialCondition which corresponds to the SweepVariable to one of the SweepValues %%%
        LocalInitialConditions(Index,ProgramParameters.sweepvar) = SweepValues(Index);
        
		%%% Create String To Use In Legend %%%
        LegendString{Index} = sprintf('%s = %.3g',Name,SweepValues(Index));

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run Simulations For All Values of SweepVariable %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor Index = 1:ProgramParameters.nsweeps

		%%% Perform ODE Calculations for the different values of SweepValue %%%
		Output{Index} = SolveODE(LocalInitialConditions(Index,:),...
						ProgramParameters, LocalModelParameters{Index});
    
end

% Testing - re-arrange structs of two elements into a single matrix
% test_Output = cell(ProgramParameters.nsweeps,1);
% for i=1:ProgramParameters.nsweeps
%     test_matrix = zeros(length(Output{i}.Y),12);
%     test_matrix(:,1:11) = Output{i}.Y;
%     test_matrix(:,12) = Output{i}.t;
%     test_Output{i} = test_matrix;
% end


    
