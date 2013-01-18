%   This script automates the compilation of the mex files in this repo.
%   "filename" is the name of the C++ file to be compiled
%   Note: assumes boost must be in the directory "boost_1_53_0_beta1" of the
%   directory this script is being run from
%   Copyright (C) 2013 Will Smith et al.
%   -----------------------------------------------------------------------
%
%   This is code is based on the work by 2013 Paul Taylor et al. on
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

close all;
clear all;

filename = 'mexodestiff1.cpp';


compileCommand = ['mex -I''' pwd '/boost_1_53_0_beta1''' ' ' filename1];
eval(compileCommand)