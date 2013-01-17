%PREPARERUN checks the validity of inputs from the GUI. Progression to the
% ODE solver SolveODE is prohibited unless the following checks are passed.
%
% Checks: min(parameter) < max(parameter) for Parameter Sweep
%			time(start) < time(end) for Simulation Parameters
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

function success = prepareRun(handles)


if ((handles.ProgramParameters.min<=handles.ProgramParameters.max) && ...
        (handles.ProgramParameters.t_start<handles.ProgramParameters.t_end))
    
	success = 1;

else 

	success = 0;

end


