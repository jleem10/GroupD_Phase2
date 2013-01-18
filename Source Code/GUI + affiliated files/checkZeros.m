%CHECKZEROS checks for non-physical results from the SolveODE and removes
% those entries. These are typically values which are below zero for
% proportions (0 to 1) or percentages (0 to 100).
%
% Example: Matrix 'eg'
% eg(any(eg<0),2) finds all rows of eg which contain an element less than
% 0. The '2' checks rows, rather than columns (which would be 1).
% eg(any(eg<0),2) = [] replaces these rows with an empty row, deleting that
% invalid set of data.
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

function result = checkZeros(result, ncells)

for Index = 1:ncells
   result{Index}.t(any(result{Index}.Y <0,2),:) = [];
   result{Index}.Y(any(result{Index}.Y <0,2),:) = [];    
end

