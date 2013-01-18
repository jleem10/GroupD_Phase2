% PLOTSWEEP contains the plotting functions (normal, semilog or loglog) for
% the sets of data produced by the parameterSweep function.
%
% It takes all the handles as inputs, and plots the correct graph in the
% appropriate set of axes.
%
% Copyright (C) 2013 Will Smith et al.
%
% -------------------------------------------------------------------------
%
%   This code is based on work by 2013 Paul Taylor et al. on
%   on reproducing the results of Hancioglu et al.'s Paper 'A Dynamical
%   Model of Human Immune Response to Influenza A Virus Infection'
%   
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
%% Size of the results matrix is changed to accommodate for new variables
%
   %%% Get results into a single matrix - easier to produce graphs %%%
    results_struct = handles.result{Index};
    results = zeros(length(results_struct.Y),15); 
    results(:,1:13) = results_struct.Y;
    results(:,15) = results_struct.t;

%
%% Scaling factor for drug interaction is incorporated
%
    if (handles.ModelParameters.t_peak == 0 || handles.ModelParameters.t_half == 0 || handles.ModelParameters.d_max == 0)
        c_1 = 0; c_2 = 0; c_3 = 0;
    else
        c_3 = 1/handles.ModelParameters.t_peak;
        c_2 = (c_3^2*handles.ModelParameters.t_half^2-4*c_3*handles.ModelParameters.t_half+1) / handles.ModelParameters.t_half;
        c_1 = handles.ModelParameters.d_max*(2*c_3+c_2);
    end;   

    
    t = results(:,15);
    results(:,14) = c_1*t./(1+c_2*t+c_3^2*t.^2);
    

%
%% Remaining portion is the same as original
%
    
    %%% Check for type of graph to be drawn %%%
    if (handles.axis1_log == [0 0])

    %%% Linear Plot %%%
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