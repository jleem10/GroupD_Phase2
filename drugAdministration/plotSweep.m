%PLOTSWEEP contains the plotting functions (normal, semilog or loglog) for
% the sets of data produced by the parameterSweep function.
%
% It takes all the handles as inputs, and plots the correct graph in the
% appropriate set of axes.
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

function plotSweep(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose Colour Scheme %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If more than one colour required, get a more interesting colour scheme
% than MATLAB's default. Otherwise, keep blue.
% Use jet colormap, and take as many colours as there were sweeps

if (handles.ProgramParameters.nsweeps == 1)
	
	colours = [0 0 1];

else

	colours = colormap(jet(handles.ProgramParameters.nsweeps));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Plot Type And Plot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Choose Whether a Normal Plot, semilog plot or loglog plot is needed %%%

for Index=1:handles.ProgramParameters.nsweeps
   
    %%% Get results into a single matrix - easier to produce graphs %%%
    results_struct = handles.result{Index};
    %%% Preallocate matrix of length of similation
    results = zeros(length(results_struct.Y),14);  %Expanded matrix to accomodate extra result
    results(:,1:11) = results_struct.Y;            %Extract Y matrix
    results(:,12) = results_struct.t;              %Extract t vector
    
    %%% Check for type of graph to be drawn %%%
    if (handles.axis1_log == [0 0])
        
        plot(results(:,handles.axis1_xaxis),results(:,handles.axis1_yaxis),...
            'Color',colours(Index,:),'LineWidth',handles.line);

    
	elseif (handles.axis1_log == [1 0])

		%%% Semi-Log Plot (X-Axis) %%%
        semilogx(results(:,handles.axis1_xaxis),results(:,handles.axis1_yaxis),...
            'Color',colours(Index,:),'LineWidth',handles.line);
    
	elseif (handles.axis1_log == [0 1])

		%%% Semi-Log Plot (Y-Axis) %%%
        semilogy(results(:,handles.axis1_xaxis),results(:,handles.axis1_yaxis),...
            'Color',colours(Index,:),'LineWidth',handles.line);
    
	else 
		
		%%% Log-Log Plot %%%
		loglog(results(:,handles.axis1_xaxis),results(:,handles.axis1_yaxis),...
			'Color',colours(Index,:),'LineWidth',handles.line);
    
	end
    
	hold on;    %%% Don't overwrite previous sets of results %%%
 
    
end

       % plot(t,results(:,13),'x') % add a reference plot for drugs
       % No lines - interpolation is poor, as
       %the solver is aiming to put most points where the coupled ODEs
       %undergo significant change - it doesn't know about the
       %concentration function!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Produce legend (if required) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (handles.ProgramParameters.nsweeps ~=1)

	legend(handles.legend_entries,'Location','EastOutside');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add labels for x and y axes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xString = getLabelNames(handles.axis1_xaxis);
yString = getLabelNames(handles.axis1_yaxis);
xlabel(xString);
ylabel(yString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Remove Hold - Allow Overwriting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold off;