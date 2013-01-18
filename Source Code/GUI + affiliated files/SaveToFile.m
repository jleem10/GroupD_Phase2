% SAVETOFILE: Function To Save Simulation Data To A Text File.
% 
% Takes in a CELL, 'Data', whose elements are STRUCTS. Each of these
% STRUCTS corresponds to the simulation data for a different value of the
% sweep parameter.
%
% i.e. Data{1} gives the simulation data for the minimum value of the sweep
% parameter.
%
% Each of the STRUCTS contained within 'Data' contains two elements, t and
% Y. Data{Index}.t is a COLUMN VECTOR containing a list of time values of
% the simulation. Data{Index}.Y is an ARRAY where each column corresponds
% to a different simulation parameter for the corresponding time value in 
% Data{Index}.t. Columns in Data{Index}.Y correspond to V(Viral), H(Healthy),
% I(Infected), M(Antigen-Presenting), F(Interferon), R(Resistant), E(Effector),
% P(Plasma), A(Antibodies), S(Antigenic Distance), D(Dead Cells).
%
% Also takes in a STRING 'FileName' which is the file to which the data will
% be saved. Should be in the form ***.txt.
%
% Copyright (C) 2013 Will Smith et al.
% 
% -------------------------------------------------------------------------
%
%   This code is based on work by 2013 Paul Taylor et al. on
%   on reproducing the results of Hancioglu et al.'s Paper 'A Dynamical
%   Model of Human Immune Response to Influenza A Virus Infection'
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

function SaveToFile( Data, FileName, SweepVariableName, SweepValues )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Open File %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen(FileName, 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write License as Header %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutputLine = 'This Influenza A Virus Infection Model Data is made available under the Open Database License: http://opendatacommons.org/licenses/odbl/1.0/.\n';
fprintf(fileID,OutputLine);
OutputLine = 'Any rights in individual contents of the database are licensed under the Database Contents License: http://opendatacommons.org/licenses/dbcl/1.0/\n\n\n';
fprintf(fileID,OutputLine);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write Variable for Reference %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutputLine = '## Variable\tDescription\n';
fprintf(fileID,OutputLine);
OutputLine = '## V\t\tViral Load\n';
fprintf(fileID,OutputLine);
OutputLine = '## H\t\tHealthy Cells\n';
fprintf(fileID,OutputLine);
OutputLine = '## I\t\tInfected Cells\n';
fprintf(fileID,OutputLine);
OutputLine = '## M\t\tActivated Antigen-Presenting Cells\n';
fprintf(fileID,OutputLine);
OutputLine = '## F\t\tInterferon\n';
fprintf(fileID,OutputLine);
OutputLine = '## R\t\tResistant Cells\n';
fprintf(fileID,OutputLine);
OutputLine = '## E\t\tEffector Cells\n';
fprintf(fileID,OutputLine);
OutputLine = '## P\t\tPlasma Cells\n';
fprintf(fileID,OutputLine);
OutputLine = '## A\t\tAntibodies\n';
fprintf(fileID,OutputLine);
OutputLine = '## S\t\tAntigenic Distance\n';
fprintf(fileID,OutputLine);
OutputLine = '## D\t\tDead Cells\n\n\n';
fprintf(fileID,OutputLine);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Print Parameter Being Varied %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutputLine = 'Sweep Variable: %s\n\n';
fprintf(fileID,OutputLine,SweepVariableName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop Over Values of Sweep Parameter and Write Data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for Index = 1:length(SweepValues)

%%% Print Current Value of Variable %%%
OutputLine = 'Sweep Variable Value: %f\n\n';
fprintf(fileID,OutputLine,SweepValues(Index));

%%% Print Table Headers %%%
OutputLine = 't\t\t\tV\t\t\tH\t\t\tI\t\t\tM\t\t\tF\t\t\tR\t\t\tE\t\t\tP\t\t\tA\t\t\tS\t\t\tD\n\n';
fprintf(fileID,OutputLine);

%%% Print Table Data %%%
OutputLine = '%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n';
fprintf(fileID,OutputLine,[Data{Index}.t,Data{Index}.Y]');

%%% Print Line Break %%%
fprintf(fileID,'\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Close File %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


