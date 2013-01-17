%   This tests the mex-file implementation of the ode solvers
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

close all;
clear all;

InitialiseModelParameters

% For the constant step size solvers this is the step size
% For the adaptive solvers this is the initial step size
options.stepSize = 0.01;

% Replace with mexode(ModelParameters, [0, 30], InitialConditions, options)
% in order to test with the constant step size solvers
tic
options.stepSize = 0.01;
[time, x] = mexodestiff(ModelParameters, [0, 30], InitialConditions, options);
toc

% Plot a phase diagram to check output
x(:,11) = 1 - x(:,2) - x(:,6) - x(:,3);
plot(x(:,11), x(:,1));

% Compare with:
% tic
% SolveODE(InitialConditions)
% toc
