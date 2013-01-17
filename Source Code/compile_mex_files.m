%   This script automates the compilation of the mex files in this project
%   "filename1" and "filename2" are the names of the C++ files to be compiled
%   Note: assumes boost must be in the directory "boost_1_52_0" of the
%   directory this script is being run from
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

close all;
clear all;

filename1 = 'mexodestiff.cpp';
filename2 = 'mexode.cpp';

compileCommand = ['mex -I''' pwd '/boost_1_52_0''' ' ' filename1];
eval(compileCommand)

compileCommand = ['mex -I''' pwd '/boost_1_52_0''' ' ' filename2];
eval(compileCommand)
