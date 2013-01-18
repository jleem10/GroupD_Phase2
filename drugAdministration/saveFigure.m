%SAVEFIGURE copies the graph on the GUI to a new figure and saves it to a
% path specified by the user. A variety of file types are available. After
% saving, the new figure is closed. This prevents the GUI from being
% included in the saved image.
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

function saveFigure(handles)

% Set size of new figure to one which fits on A4
% Make a wider figure if legens is to be printed too
if (handles.ProgramParameters.nsweeps == 1)
    handles.save = figure('Units','centimeters','Position',[5 5 16 14]);
else handles.save = figure('Units','centimeters','Position',[5 5 20 14]);
end

set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPositionMode','auto');

% Re-plot graph
if (handles.success == 1)
plotSweep(handles);
else return
end

% Change x and y axis limits
xLimits = [str2double(get(handles.edit_plot_minx,'String')),...
    str2double(get(handles.edit_plot_maxx,'String'))];
yLimits = [str2double(get(handles.edit_plot_miny,'String')),...
    str2double(get(handles.edit_plot_maxy,'String'))];
xlim(xLimits);
ylim(yLimits);

str2double(get(handles.edit_plot_minx,'String'));

% Choose filename, path and file type with new input box
[filename, pathname, filterindex] = uiputfile({'*.pdf';'*.png';'*.jpg';},...
    'Save graph as...');
    
% Can exit dialog box without saving
if isequal(filename,0) || isequal(pathname,0)
    return
end

% Save figure with chosen attributes, then close new figure
saveas(gcf,fullfile(pathname,filename))
close(handles.save)

% Return to main GUI figure
figure(handles.figure1);