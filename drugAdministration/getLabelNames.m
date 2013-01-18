%GETLABELNAMES chooses the axis label names from a list of possible
% options, based on the current selection of the x and y popupmenus in the
% plotting section.
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

function name = getLabelNames(number)

switch number
    case 1
        name = 'V';
    case 2
        name = 'H';
    case 3
        name = 'I';
    case 4
        name = 'M';
    case 5
        name = 'F';
    case 6
        name = 'R';
    case 7
        name = 'E';
    case 8
        name = 'P';
    case 9
        name = 'A';
    case 10
        name = 'S';
    case 11
        name = 'D';
    case 12
        name = 't / days';
end

